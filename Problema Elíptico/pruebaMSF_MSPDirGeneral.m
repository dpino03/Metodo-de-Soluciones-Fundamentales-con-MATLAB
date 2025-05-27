% -------------------------------------------------------------------------
% Veamos un ejemplo de implementación del código MSFMSPDirGeneral para
% valores concretos.
% -------------------------------------------------------------------------

% Se presentan varias opciones de dominios para probar.

% % Dominio con forma de estrella
% xc = [0;0];
% k = 5; % Número de picos
% x1_eq = @(theta) (1 + 0.3*cos(k*theta)).*cos(theta);
% x2_eq = @(theta) (1 + 0.3*cos(k*theta)).*sin(theta);

% Dominio con forma de rosa polar:

xc = [0; 0];           % Centroide (origen)
n_pet = 7;             % Número de pétalos 
r0 = 1;                % Radio base
amp = 0.5;             % Amplitud de oscilación de pétalos
x1_eq = @(theta) (r0 + amp*cos(n_pet*theta)) .* cos(theta);
x2_eq = @(theta) (r0 + amp*cos(n_pet*theta)) .* sin(theta);

dist = .4;      % Distancia para el mallado interior 
Nb = 100;       % Número de puntos sobre Gamma1
lambda0 = 1.5;  % Parámetro para la colocación de las fuentes 
lambda1 = 10;   % Parámetro para definir la función f
 
u_exacta = @(x) exp(x(2))*sin(x(1));                 % Solución exacta
a = @(x) sqrt(abs(x(1)+x(2)));                       % Función a(x)
h = @(x) sqrt(abs(x(1)+x(2)))*exp(x(2))*sin(x(1));   % Función h(x)
g = @(x) exp(x(2))*sin(x(1));                        % Condición Dirichlet

% Resolución:

[u] = MSF_MSPDirGeneral(x1_eq,x2_eq,xc,dist,Nb,lambda0,lambda1,a,h,g,u_exacta);