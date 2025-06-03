function [u] = MSF_MSPDirGeneral(x1_eq,x2_eq,xc,dist,Nb,lambda0,lambda1,a,h,g,u_exacta)
% -------------------------------------------------------------------------
% Resolución del problema: 
%
% -Lap(u) + a*u = h en Omega
%  u = g            sobre fr(Omega)
%
% Entradas:
% x1_eq, x2_eq: ecuaciones paramétricas de la frontera.
% xc: centroide del dominio.
% dist: resolución del mallado interno.
% Nb: número de puntos frontera/fuentes.
% lambda0: factor para colocar las fuentes fuera del dominio.
% lambda1: parámetro en f y F (Lap(F) = f).
% a, h: funciones del problema.
% g: condición de Dirichlet.
% u_exacta: solución exacta (para comparar).
%
% Salida:
%
% Aproximación de la solución u(x) por el MSF-MSP y representaciones
% gráficas de la solución del método, la exacta y el error.
% -------------------------------------------------------------------------

% Discretización de la frontera.

theta = linspace(0,2*pi,Nb+1); theta(end) = [];
x1_b = x1_eq(theta);
x2_b = x2_eq(theta);
p_frontera = [x1_b;x2_b];

%  Puntos fuente fuera del dominio.

p_fuentes = p_frontera + lambda0 * (p_frontera - xc);

% Construcción del mallado interior.

x1_min = min(x1_b); x1_max = max(x1_b);
x2_min = min(x2_b); x2_max = max(x2_b);
[x1g,x2g] = meshgrid(x1_min:dist:x1_max,x2_min:dist:x2_max);
in = inpolygon(x1g,x2g,x1_b,x2_b);
x1_in = x1g(in); x2_in = x2g(in);
p_interior = [x1_in'; x2_in'];

Nf = size(p_interior,2);

% Definición  de las funciones necesarias.

G = @(r) -1/(2*pi)*log(max(r,1e-10)); 
f = @(r) (r <= lambda1).*(1 - max(r,1e-10)/lambda1).^2;
F = @(r) (r <= lambda1) .* (r.^4/(16*lambda1^2) - (2*r.^3)/(9*lambda1) + r.^2/4) + ...
         (r > lambda1) .* (13*lambda1^2 + lambda1^2/12*log(max(r/lambda1,1e-10)));

% Construcción de las matrices y el término independiente para el sistema.

C = zeros(Nf,Nf);
for i = 1:Nf
    ri = p_interior(:,i);
    for j = 1:Nf
        rj = p_interior(:,j);
        C(i,j) = -f(norm(ri - rj)) + a(ri)*F(norm(ri - rj));
    end
end

D = zeros(Nf,Nb);
for i = 1:Nf
    ri = p_interior(:,i);
    for j = 1:Nb
        sj = p_fuentes(:,j);
        D(i,j) = a(ri) * G(norm(ri - sj));
    end
end

E = zeros(Nb,Nf);
for i = 1:Nb
    xi = p_frontera(:,i);
    for j = 1:Nf
        rj = p_interior(:,j);
        E(i,j) = F(norm(xi - rj));
    end
end

L = zeros(Nb,Nb);
for i = 1:Nb
    xi = p_frontera(:,i);
    for j = 1:Nb
        sj = p_fuentes(:,j);
        L(i,j) = G(norm(xi - sj));
    end
end

A = [C,D;E,L];
lambda = 1e-6;
A = A + lambda*eye(size(A));

b1 = zeros(1,Nf);
for i = 1:Nf
    ri = p_interior(:,i);
    b1(i) = h(ri);
end

b2 = zeros(1,Nb);
for i = 1:Nb
    x_i = p_frontera(:,i);
    b2(i) = g(x_i);
end

b = [b1,b2];
b = b';

% Resolución del sistema.

coefs = A\b;
beta = coefs(1:Nf);
alpha = coefs(Nf+1:end);

% Definición de la solución u(x).

u = @(x) dot(beta,arrayfun(@(j) F(norm(x - p_interior(:,j))), 1:Nf)) + ...
         dot(alpha,arrayfun(@(k)G(norm(x - p_fuentes(:,k))), 1:Nb));

% Evaluación en una malla para la visualización.

Nmalla = 400;
[X,Y] = meshgrid(linspace(x1_min,x1_max,Nmalla),linspace(x2_min,x2_max,Nmalla));
in = inpolygon(X,Y,x1_b,x2_b);

U = nan(size(X));
idx = find(in);
[i_idx,j_idx] = ind2sub(size(in),idx);
U(idx) = arrayfun(@(i,j) u([X(i,j); Y(i,j)]),i_idx,j_idx);

% Evaluación en la malla de la solución exacta.

Ue = nan(size(X));
Ue(idx) = arrayfun(@(i,j) u_exacta([X(i,j);Y(i,j)]),i_idx,j_idx);

% Gráficas de los resultados.

figure(1);
surf(X,Y,U,'EdgeColor','none'); view(3); colorbar;
xlabel('x_1');
ylabel('x_2');
title('Solución del MSF');

figure(2);
surf(X,Y,Ue,'EdgeColor','none'); view(3); colorbar;
xlabel('x_1');
ylabel('x_2');
title('Solución exacta');

figure(3);
surf(X,Y,abs(U - Ue),'EdgeColor','none'); view(3); colorbar;
xlabel('x_1');
ylabel('x_2');
title('Error absoluto');

% Impresión del error absoluto máximo.

disp(['Error absoluto máximo: ',num2str(max(abs(U(:)-Ue(:))))]);

% Gráficar de la distribución de los puntos y las fuentes en el dominio.

figure(4);
hold on;
scatter(p_interior(1,:), p_interior(2,:), 10, 'b', 'filled');
scatter(p_frontera(1,:), p_frontera(2,:), 30, 'r', 'filled');
scatter(p_fuentes(1,:), p_fuentes(2,:), 50, 'm', 'x', 'LineWidth', 1.5);
plot(x1_eq(theta), x2_eq(theta), 'k-', 'LineWidth', 1.5);
legend('Interior','Frontera','Fuentes','Frontera real');
axis equal; grid on; title('Distribución de puntos y fuentes');
xlabel('x_1'); ylabel('x_2');

end
