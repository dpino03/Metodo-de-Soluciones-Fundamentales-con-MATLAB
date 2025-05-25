function [u] = MSFLaplaceGeneral(x_eq,y_eq,xc,N,lambda,g,u_exacta)
% -------------------------------------------------------------------------
% Resolución del problema de Laplace en un dominio general 2D definido por
% ecuaciones paramétricas.
%
% Lap(u) = 0 en Omega
% u = g      sobre fr(Omega)
%
% Entradas:
% 
% x_eq = Función handle para x(theta), theta en [0,2pi].
% y_eq = Función handle para y(theta), theta en [0,2pi].
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
puntos_frontera = [x_eq(theta);y_eq(theta)];

% Colocamos las fuentes en la frontera virtual.

puntos_fuentes = puntos_frontera + lambda*(puntos_frontera-xc);

% Construimos el sistema.

G = @(x,y) -1/(2*pi)*log(sqrt(x.^2+y.^2+1e-10)); % Solución fundamental

A = zeros(N, N);

for i = 1:N
    r_x = puntos_frontera(1,i) - puntos_fuentes(1,:);
    r_y = puntos_frontera(2,i) - puntos_fuentes(2,:);
    A(i,:) = G(r_x,r_y);
end

b = g(puntos_frontera(1,:), puntos_frontera(2,:))';

% Resolvemos el sistema usando un método de regularización.

At = transpose(A);
epsilon = 1e-8; % Parámetro de regularización
B = At*A + epsilon*eye(size(A));
d = At*b;
coefs = B\d;

% Definimos la solución del método.

u = @(x,y) sum(coefs'.* G(x-puntos_fuentes(1,:),y-puntos_fuentes(2,:)),2);

% Finalmente, vamos a crear una malla y representar gráficamente la
% solución.

Nmalla = 1000;

theta_eval = linspace(0,2*pi,Nmalla);
x_front = x_eq(theta_eval);
y_front = y_eq(theta_eval);

x_min = min(x_front); x_max = max(x_front);
y_min = min(y_front); y_max = max(y_front);

% Creamos una malla adaptada al dominio.

[X,Y] = meshgrid(linspace(x_min, x_max, Nmalla), ...
                  linspace(y_min, y_max, Nmalla));

% Determinamos  qué puntos están dentro del dominio.

in = inpolygon(X(:),Y(:),x_front,y_front);
U = nan(size(X));
U(in) = u(X(in),Y(in));

% Procedemos con la visualización de la solución del método.

figure(1);
surf(X, Y, U,'EdgeColor','none');
hold on;
xlabel('x'); ylabel('y'); zlabel('u(x,y)');
title('Solución del MSF');
colorbar;
view(2);
axis equal;

% Comparamos con la solución exacta y calculamos el error.

U_exact = nan(size(X));
U_exact(in) = u_exacta(X(in),Y(in));
    
figure(2);
surf(X,Y,U_exact,'EdgeColor','none');
title('Solución exacta');
view(2); 
colorbar;
axis equal;
    
figure(3);
error = abs(U-U_exact);
surf(X,Y,error,'EdgeColor','none');
title('Error absoluto');
view(2); 
colorbar; 
axis equal;


% -------------------------------------------------------------------------
% Graficar puntos del mallado, frontera y fuentes.
% -------------------------------------------------------------------------

figure(4);
hold on;

% 2. Puntos en la frontera (Dirichlet, rojo)
scatter(puntos_frontera(1,:), puntos_frontera(2,:), 30, 'r', 'filled');

% 4. Fuentes fuera de la bola (magenta)
scatter(puntos_fuentes(1,:), puntos_fuentes(2,:), 50, 'm', 'x', 'LineWidth', 1.5);

% Dibujar la frontera del dominio
theta_rep = linspace(0,2*pi,100);
plot(x_eq(theta_rep),y_eq(theta_rep), 'k-', 'LineWidth', 1.5);

% Configuración de la gráfica
xlabel('x_1', 'Color', 'k', 'FontWeight', 'bold');
ylabel('x_2', 'Color', 'k', 'FontWeight', 'bold');
title('Distribución de puntos y fuentes', 'Color', 'k', 'FontWeight', 'bold');
legend('Puntos sobre la frontera', 'Fuentes', 'Location', 'best');
grid on;
axis equal;
hold off;

saveas(figure(1),'SolMetLaplaceGeneral.png');
saveas(figure(2),'SolExactaLaplaceGeneral.png');
saveas(figure(3),'ErrorLaplaceGeneral.png');
saveas(figure(4),'puntosGeneral.png');

end

