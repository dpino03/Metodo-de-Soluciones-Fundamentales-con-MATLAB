% -------------------------------------------------------------------------
% Veamos un ejemplo de implementación del código MSFCalor para
% valores concretos.
% -------------------------------------------------------------------------

R = 1;          % Radio de la bola
xc = [0,0];     % Centro de la bola
T = 1;          % Tiempo final
M = 50;         % Número de subintervalos temporales
dist = 0.2;     % Distancia entre puntos del mallado regular
Nb = 50;        % Número de puntos en la frontera
lambda0 = 1.5;  % Factor para ubicar fuentes fuera de la bola
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

[u] = MSFCalor(R,xc,T,M,u0,dist,Nb,lambda0,lambda1,h,g,u_exacta);



