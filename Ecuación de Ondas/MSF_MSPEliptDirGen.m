function [u] = MSF_MSPEliptDirGen(x1_eq,x2_eq,xc,dist,Nb,lambda0,lambda1,a,h,g)
% -------------------------------------------------------------------------
% Consideramos un problema de la forma:
%
% -Lap(u) + a*u = h en Omega
% u = g             sobre fr(Omega)
%
% Entradas:
%
% x1_eq, x2_eq: ecuaciones paramétricas de la frontera.
% xc = Centro de la bola (vector 2x1).
% dist = Distancia entre los puntos del mallado regular interior para el
% cálculo de la solución particular u_P.
% Nb = Número de puntos de la frontera para aplicar el MSF
% lambda0 = Factor para ubicar las fuentes fuera de la bola.
% lambda1 = Constante que aparece en la definición de las funciones
% radiales f y F, relacionadas por Lap(F) = f.
% a = Función que describe la EDP.
% h = Función que describe la EDP.
% g = Función que describe las condiciones de Dirichlet en la frontera.
% u_exacta = Solución exacta para realizar la comparación
%
% Salida:
%
% Aproximación de la solución u(x) por el MSF-MSP.
% -------------------------------------------------------------------------

% Discretización de la frontera.

theta = linspace(0,2*pi,Nb+1); theta(end) = [];
x1_b = x1_eq(theta);
x2_b = x2_eq(theta);
p_frontera = [x1_b;x2_b];

% Puntos fuente fuera del dominio.

p_fuentes = p_frontera + lambda0 * (p_frontera - xc);

% Construcción del mallado interior.

x1_min = min(x1_b); x1_max = max(x1_b);
x2_min = min(x2_b); x2_max = max(x2_b);
[x1g,x2g] = meshgrid(x1_min:dist:x1_max,x2_min:dist:x2_max);
in = inpolygon(x1g,x2g,x1_b,x2_b);
x1_in = x1g(in); x2_in = x2g(in);
p_interior = [x1_in'; x2_in'];

Nf = size(p_interior,2);

% Definición de funciones necesarias.

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