% -------------------------------------------------------------------------
% Veamos un ejemplo de implementación del código MSFLaplaceBola para
% valores concretos.
% -------------------------------------------------------------------------

R = 1;                               % Radio de la bola
xc = [0;0];                          % Centro de la bola
N = 100;                             % Número de puntos en la frontera
lambda = 1.75;                       % No se puede tomar muy pequeño
g = @(x) exp(x(1))*cos(x(2));        % Condición de contorno
u_exacta = @(x) exp(x(1))*cos(x(2)); % Solución exacta

% Nota: Si cogemos demasiados puntos (N grande), la matriz estará mal 
% condicionada, y el coste computacional será bastante elevado.

% Nota: Si cogemos lambda muy pequeño, la solución se distorsionará, puesto
% que la solución phi(x,s) tiene una singularidad en x = s.

% Llamamos a la función y comprobamos un valor de prueba.

[u] = MSFLaplaceBola(R,xc,N,lambda,g,u_exacta);

prueba_valor_concreto = u([pi;pi]);

disp(['El valor de prueba es ',num2str(prueba_valor_concreto),'.'])

% Nota: En el caso de que la función tenga varios puntos donde se anule, el error
% relativo no es representativo, ya que hay varios puntos del dominio con 
% valores de la solución exacta cercanos a 0. Sin embargo, vemos con el 
% error absoluto que la aproximación es muy buena.