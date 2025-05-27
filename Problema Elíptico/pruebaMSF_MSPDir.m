% -------------------------------------------------------------------------
% Veamos un ejemplo de implementación del código MSFMSPDir para
% valores concretos.
% -------------------------------------------------------------------------

R = 1;          % Radio de la bola
xc = [0;0];     % Centro de la bola
dist = .1;      % Distancia para el mallado interior 
Nb = 100;       % Número de puntos sobre Gamma1
lambda0 = 1.5;  % Parámetro para la colocación de las fuentes 
lambda1 = 10;   % Parámetro para definir la función f
 
u_exacta = @(x) exp(x(2))*sin(x(1));                 % Solución exacta
a = @(x) sqrt(abs(x(1)+x(2)));                       % Función a(x)
h = @(x) sqrt(abs(x(1)+x(2)))*exp(x(2))*sin(x(1));   % Función h(x)
g = @(x) exp(x(2))*sin(x(1));                        % Condición Dirichlet

% Llamamos a la función y comprobamos un valor de prueba.

[u] = MSF_MSPDir(R,xc,dist,Nb,lambda0,lambda1,a,h,g,u_exacta);

valor_prueba = u([1;1]);

disp(['El valor de prueba es ',num2str(valor_prueba)]);

% Nota: Hay que ajustar los parámetros lambda0 y lambda1 para obtener una
% buena aproximación.

% Nota: Para el problema con condiciones sólo de  tipo Dirichlet, es mejor 
% que lambda1 no tome valores demasiado elevados, en este caso, tomamos
% lambda1 = 10.


