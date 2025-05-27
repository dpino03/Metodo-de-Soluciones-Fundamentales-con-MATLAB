% -------------------------------------------------------------------------
% Veamos un ejemplo de implementación del código MSFLaplaceGeneral para
% valores concretos.
% -------------------------------------------------------------------------

% Dominio con forma de estrella:

xc = [0;0];
k = 5; % Número de picos
x1_eq = @(theta) (1 + 0.3*cos(k*theta)).*cos(theta);
x2_eq = @(theta) (1 + 0.3*cos(k*theta)).*sin(theta);

g = @(x) exp(x(1)).*cos(x(2));        % Condición de contorno
u_exacta = @(x) exp(x(1)).*cos(x(2));

% Parámetros:

N = 100;       % Puntos frontera
lambda = 1.5;  % Distancia para fuentes

% Resolución:

[u] = MSFLaplaceGeneral(x1_eq,x2_eq,xc,N,lambda,g,u_exacta);