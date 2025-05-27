% -------------------------------------------------------------------------
% Veamos un ejemplo de implementación del código MSFOndasGeneral para
% valores concretos.
% -------------------------------------------------------------------------

k = 3; 
x1_eq = @(theta) (1 + 0.3*cos(k*theta)).*cos(theta);
x2_eq = @(theta) (1 + 0.3*cos(k*theta)).*sin(theta);
xc = [0;0];     % Centroide
T = 1;          % Tiempo final
M = 50;         % Número de subintervalos temporales
dist = 0.2;     % Distancia entre puntos del mallado regular
Nb = 50;        % Número de puntos en la frontera
lambda0 = 1.75; % Factor para ubicar fuentes fuera de la bola
lambda1 = 10;   % Constante para las funciones radiales

% Función h(x,t) (término fuente).

h = @(x,t) (2*pi^2/100)*(t+1)*cos(pi*x(1)/10)*cos(pi*x(2)/10);

% Función g(x,t) (condición de Dirichlet en la frontera).

g = @(x,t) cos(pi*x(1)/10)*cos(pi*x(2)/10)*(t+1);

% Condición inicial para u, u0(x).

u0 = @(x) cos(pi*x(1)/10)*cos(pi*x(2)/10);

% Condición inicial para u_t, v0(x).

v0 = @(x) cos(pi*x(1)/10)*cos(pi*x(2)/10);

% Solución exacta u_exacta(x,t).

u_exacta = @(x,t) cos(pi*x(1)/10)*cos(pi*x(2)/10)*(t+1);

% Implementación.

[u] = MSFOndasGeneral(x1_eq,x2_eq,xc,T,M,u0,v0,dist,Nb,lambda0,lambda1,h,g,u_exacta);