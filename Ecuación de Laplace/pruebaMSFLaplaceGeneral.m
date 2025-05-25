% Dominio con forma de estrella
xc = [0;0];
k = 5; % Número de picos
x_eq = @(theta) (1 + 0.3*cos(k*theta)).*cos(theta);
y_eq = @(theta) (1 + 0.3*cos(k*theta)).*sin(theta);

g = @(x,y) exp(x).*cos(y);        % Condición de contorno
u_exacta = @(x,y) exp(x).*cos(y);

% Parámetros:

N = 100;       % Puntos frontera
lambda = 1.75; % Distancia para fuentes

% Resolución:

[u] = MSFLaplaceGeneral(x_eq,y_eq,xc,N,lambda,g,u_exacta);