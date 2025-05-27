function [u] = MSFLaplaceGeneral(x1_eq,x2_eq,xc,N,lambda,g,u_exacta)
% -------------------------------------------------------------------------
% Resolución del problema de Laplace en un dominio general 2D definido por
% ecuaciones paramétricas.
%
% Lap(u) = 0 en Omega
% u = g      sobre fr(Omega)
%
% Entradas:
% 
% x1_eq = Función handle para x1(theta), theta en [0,2pi].
% x2_eq = Función handle para x2(theta), theta en [0,2pi].
% xc = Centro geométrico o centroide.
% N = Número de puntos de frontera.
% lambda = Factor de distancia para fuentes.
% g = Función de condición de Dirichlet en la frontera: g(x,y).
% u_exacta = Solución exacta para la comparación.
%
% Salidas:
% u = Función la solución aproximada
% Representación gráfica de la solución exacta, la del método, y el error.
% -------------------------------------------------------------------------

% Comenzamos discretizando la frontera.

theta = linspace(0,2*pi,N+1);
theta(end) = [];
puntos_frontera = [x1_eq(theta);x2_eq(theta)];

% Colocamos las fuentes en la frontera virtual.

puntos_fuentes = puntos_frontera + lambda*(puntos_frontera-xc);

% Construimos el sistema.

G = @(r) -1/(2*pi)*log(max(norm(r),1e-10));% Solución fundamental

A = zeros(N,N);

for i = 1:N
    x_i = puntos_frontera(:,i);    % Punto en la frontera x_i
    for j = 1:N
        s_j = puntos_fuentes(:,j); % Fuente s_j
        A(i,j) = G(x_i-s_j);       % Evaluación de la solución fundamental
    end
end

b = zeros(N,1);

for i = 1:N
    b(i) = g(puntos_frontera(:,i)); % Condiciones de Dirichlet
end

% Resolvemos el sistema usando un método de regularización.

At = transpose(A);
epsilon = 1e-8; % Parámetro de regularización
B = At*A + epsilon*eye(size(A));
d = At*b;
coefs = B\d;

% Definimos la solución del método.

u = @(x) dot(coefs,arrayfun(@(j) G(x-puntos_fuentes(:,j)),1:N));

% Finalmente, vamos a crear una malla y representar gráficamente la
% solución.

Nmalla = 500;

theta_eval = linspace(0,2*pi,Nmalla);
x1_front = x1_eq(theta_eval);
x2_front = x2_eq(theta_eval);

x1_min = min(x1_front); x1_max = max(x1_front);
x2_min = min(x2_front); x2_max = max(x2_front);

% Creamos una malla adaptada al dominio.

[X,Y] = meshgrid(linspace(x1_min,x1_max,Nmalla), ...
                  linspace(x2_min,x2_max,Nmalla));

% Determinamos  qué puntos están dentro del dominio.

in = inpolygon(X(:),Y(:),x1_front,x2_front);
U = nan(size(X));
idx = find(in);
U(idx) = arrayfun(@(k) u([X(idx(k)); Y(idx(k))]), 1:length(idx));

% Procedemos con la visualización de la solución del método.

figure(1);
surf(X,Y,U,'EdgeColor','none');
hold on;
xlabel('x_1'); ylabel('x_2'); zlabel('u(x_1,x_2)');
title('Solución del MSF');
colorbar;
view(3);
axis equal;

% Comparamos con la solución exacta y calculamos el error.

Ue = nan(size(X));
Ue(idx) = arrayfun(@(k) u_exacta([X(idx(k)); Y(idx(k))]), 1:length(idx));
    
figure(2);
surf(X,Y,Ue,'EdgeColor','none');
title('Solución exacta');
view(3); 
colorbar;
axis equal;
    
figure(3);
error = abs(U-Ue);
surf(X,Y,error,'EdgeColor','none');
title('Error absoluto');
view(3); 
colorbar; 
axis equal;


% -------------------------------------------------------------------------
% Graficar puntos del mallado, frontera y fuentes.
% -------------------------------------------------------------------------

figure(4);
hold on;

% Puntos en la frontera (Dirichlet, rojo).
scatter(puntos_frontera(1,:), puntos_frontera(2,:), 30, 'r', 'filled');

% Fuentes fuera de la bola (magenta).
scatter(puntos_fuentes(1,:), puntos_fuentes(2,:), 50, 'm', 'x', 'LineWidth', 1.5);

% Frontera del dominio.
theta_rep = linspace(0,2*pi,100);
plot(x1_eq(theta_rep),x2_eq(theta_rep), 'k-', 'LineWidth', 1.5);

% Configuración de la gráfica.
xlabel('x_1', 'Color', 'k', 'FontWeight', 'bold');
ylabel('x_2', 'Color', 'k', 'FontWeight', 'bold');
title('Distribución de puntos y fuentes', 'Color', 'k', 'FontWeight', 'bold');
legend('Puntos sobre la frontera', 'Fuentes', 'Location', 'best');
grid on;
axis equal;
hold off;

end

