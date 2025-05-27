% -------------------------------------------------------------------------
% Veamos un ejemplo de implementación del código PbInversoD para
% valores concretos.
% -------------------------------------------------------------------------

R = 2;          % Radio de la bola
r = 1;          % Radio pequeño.
xc = [0;0];     % Centro de la bola
dist = .4;      % Distancia para el mallado interior 
Nb = 50;        % Número de puntos sobre Gamma1
lambda0 = 1.75; % Parámetro para la colocación de las fuentes 
lambda1 = 10;   % Parámetro para definir la función f
 
A = @(x) x(1)+x(2);                       
h = @(x) (abs(x(1)+x(2)))*exp(x(2))*sin(x(1)); 
g = @(x) exp(x(2))*sin(x(1)); 

pert = 1.e-2;
a0 = [0,0,0];

% Resolución:

[a_opt] = PbInversoD(R,r,xc,dist,Nb,lambda0,lambda1,A,h,g,pert,a0);

fprintf('El valor del polinomio óptimo es %.2f + %.2fx(1) + %.2fx(2).\n', ...
         a_opt(1),a_opt(2),a_opt(3));


% Ejemplo: 

% Con A = x(1)+x(2), pert = 1.e-4, a0 = [0,0,0] convergencia en 16
% iteraciones a -0.00 + 1.00*x(1) + 1.00*x(2). V.o: 1.755547e-10. 

% Con A = x(1)+x(2), pert = 1.e-3, a0 = [0,0,0] convergencia en 16
% iteraciones a -0.00 + 1.00*x(1) + 1.00*x(2). V.o: 1.471382e-08. 

% Con A = x(1)+x(2), pert = 1.e-2, a0 = [0,0,0] convergencia en 16
% iteraciones a 0.02 + 1.01*x(1) + 1.01*x(2). V.o: 1.229420e-06. 




