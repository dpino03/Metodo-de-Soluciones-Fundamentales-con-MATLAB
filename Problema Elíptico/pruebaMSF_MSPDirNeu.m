% -------------------------------------------------------------------------
% Veamos un ejemplo de implementación del código MSFMSPDirNeu para
% valores concretos.
% -------------------------------------------------------------------------

R = 1;          % Radio de la bola
xc = [0;0];     % Centro de la bola
dist = .2;      % Distancia para el mallado interior 
Nb1 = 50;       % Número de puntos sobre Gamma1
Nb2 = 50;       % Número de puntos sobre Gamma2
lambda0 = 1.5;  % Parámetro para la colocación de las fuentes 
lambda1 = 1.e8; % Parámetro para definir la función f
 
u_exacta = @(x) exp(x(2))*sin(x(1));                 % Solución exacta
a = @(x) sqrt(abs(x(1)+x(2)));                       % Función a(x)
h = @(x) sqrt(abs(x(1)+x(2)))*exp(x(2))*sin(x(1));   % Función h(x)
g = @(x) exp(x(2))*sin(x(1));                        % Condición Dirichlet
xi = @(x) (exp(x(2))*cos(x(1))*(x(1)-xc(1))/R + ...
          exp(x(2))*sin(x(1))*(x(2)-xc(2)))/R;       % Condición Neumann

% Llamamos a la función y comprobamos un valor de prueba.

[u] = MSF_MSPDirNeu(R,xc,dist,Nb1,Nb2,lambda0,lambda1,a,h,g,xi,u_exacta);

valor_prueba = u([1;1]);

disp(['El valor de prueba es ',num2str(valor_prueba)]);

% Nota: Hay que ajustar los parámetros lambda0 y lambda1 para obtener una
% buena aproximación.

% Nota: Tomando lambda1 arbitrariamente grande se consigue una mejor aproximación
% de la solución.

