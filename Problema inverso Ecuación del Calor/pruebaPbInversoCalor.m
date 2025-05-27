R = 1;          % Radio de la bola
xc = [0,0];     % Centro de la bola
T = 1;          % Tiempo final
M = 20;         % Número de subintervalos temporales
dist = 0.4;     % Distancia entre puntos del mallado regular
Nb = 50;        % Número de puntos en la frontera
lambda0 = 1.75; % Factor para ubicar fuentes fuera de la bola
lambda1 = 10;   % Constante para las funciones radiales

H = @(x,t) 10;
g = @(x,t) exp(-t)*sin(pi*x(1)).*sin(pi*x(2)); 
u0 = @(x) sin(pi*x(1))*sin(pi*x(2));   

pert = 1.e-1;
h0 = 0;

% Resolución:

[a_opt] = PbInversoCalor(R,xc,T,M,u0,dist,Nb,lambda0,lambda1,H,g,pert,h0);

fprintf('El valor de h óptimo es la constante %.2f\n',a_opt);


% Ejemplo:

% Para H = 10 y pert = 1.e-4 -> 3 iter, devuelve 10.00. V.o.: 3.109938e-06.

% Para H = 10 y pert = 1.e-3 -> 4 iter, devuelve 10.00. V.o.: 2.582204e-04.

% Para H = 10 y pert = 1.e-2 -> 4 iter, devuelve 10.00. V.o.: 2.809573e-02.

% Para H = 10 y pert = 1.e-1 -> 4 iter, devuelve 9.96. V.o.: 2.772585e+00.



