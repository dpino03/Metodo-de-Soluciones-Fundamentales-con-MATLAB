% -------------------------------------------------------------------------
% Veamos un ejemplo de implementación del código MSFCalorGeneral para
% valores concretos.
% -------------------------------------------------------------------------

k = 3; 
x1_eq = @(theta) (1 + 0.3*cos(k*theta)).*cos(theta);
x2_eq = @(theta) (1 + 0.3*cos(k*theta)).*sin(theta);
xc = [0;0];     % Centro geométrico
T = 1;          % Tiempo final
M = 50;         % Número de subintervalos temporales
dist = 0.1;     % Distancia entre puntos del mallado regular
Nb = 100;       % Número de puntos en la frontera
lambda0 = 1.75; % Factor para ubicar fuentes fuera de la bola
lambda1 = 10;   % Constante para las funciones radiales

% Función f(x,t) (término fuente).

h = @(x,t) (-1+2*pi^2)*exp(-t)*sin(pi*x(1))*sin(pi*x(2));

% Función g(x,t) (condición de Dirichlet en la frontera).

g = @(x,t) exp(-t)*sin(pi*x(1)).*sin(pi*x(2)); 

% Condición inicial u0(x).

u0 = @(x) sin(pi*x(1))*sin(pi*x(2)); % Condición inicial en t = 0

% Solución exacta u_exacta(x,t).

u_exacta = @(x,t) exp(-t)*sin(pi*x(1))*sin(pi*x(2)); % Solución exacta conocida

% Implementación.

[u] = MSFCalorGeneral(x1_eq,x2_eq,xc,T,M,u0,dist,Nb,lambda0,lambda1,h,g,u_exacta);