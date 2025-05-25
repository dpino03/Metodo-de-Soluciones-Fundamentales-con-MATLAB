% -------------------------------------------------------------------------
% Veamos un ejemplo de implementación del código MSFLaplaceBola para
% valores concretos.
% -------------------------------------------------------------------------

l1 = 5;                              % Longitud lados horizontales
l2 = 5;                              % Longitud lados verticales
xc = [0;0];                          % Centro de la bola
N = 100;                             % Número de puntos en la frontera
lambda = 1.75;                       % No se puede tomar muy pequeño
g = @(x) exp(x(1))*cos(x(2));        % Condición de contorno
u_exacta = @(x) exp(x(1))*cos(x(2)); % Solución exacta

% Llamamos a la función y comprobamos un valor de prueba.

[u] = MSFLaplaceRectangulo(xc,l1,l2,N,lambda,g,u_exacta);

prueba_valor_concreto = u([.2;.3]);

disp(['El valor de prueba es ',num2str(prueba_valor_concreto),'.'])
