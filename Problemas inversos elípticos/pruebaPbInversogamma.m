% -------------------------------------------------------------------------
% Veamos un ejemplo de implementación del código PbInversogamma para
% valores concretos.
% -------------------------------------------------------------------------

R = 1;          % Radio de la bola
xc = [1;2];     % Centro de la bola
dist = .2;      % Distancia para el mallado interior 
Nb = 50;        % Número de puntos sobre Gamma1
lambda0 = 1.75; % Parámetro para la colocación de las fuentes 
lambda1 = 10;   % Parámetro para definir la función f
 
A = @(x) x(1)+x(2);                       
h = @(x) sqrt(abs(x(1)+x(2)))*exp(x(2))*sin(x(1)); 
g = @(x) exp(x(2))*sin(x(1)); 

pert = 1.e-3;
a0 = [0,0,0];

% Resolución:

[a_opt] = PbInversogamma(R,xc,dist,Nb,lambda0,lambda1,A,h,g,pert,a0);

fprintf('El valor del polinomio óptimo es %.2f + %.2fx(1) + %.2fx(2).\n', ...
         a_opt(1),a_opt(2),a_opt(3));

% Ejemplos:

% Con A = 1+2*x(1)+4*x(2), pert = 1.e-4, a0 = [0,0,0], convergencia en 18
% iteraciones a 1.00 + 2.00*x(1) + 4.00*x(2). V.o: 3.644543e-06.

% Con A = 1+2*x(1)+4*x(2), pert = 1.e-3, a0 = [0,0,0] convergencia en 17
% iteraciones a 1.11 + 2.00x(1) + 3.95x(2). V.o:  2.965464e-04. 

% Con A = 1+2*x(1)+4*x(2), pert = 1.e-2, a0 = [0,0,0] convergencia en 18
% iteraciones a 0.86 + 1.97x(1) + 4.05x(2). V.o: 1.644015e-01 . 

%-----------

% Con A = x(1)+x(2), pert = 1.e-4, a0 = [0,0,0] convergencia en 16
% iter a 0.00 + 1.00x(1) + 1.00x(2). V.o: 2.611114e-07. 

% Con A = x(1)+x(2), pert = 1.e-3, a0 = [0,0,0] convergencia en 17
% iteraciones a -0.01 + 1.00x(1) + 1.00x(2). V.o: 1.203192e-04.

% Con A = x(1)+x(2), pert = 1.e-2, a0 = [0,0,0] convergencia en 15 
% iteraciones a 0.11 + 1.02x(1) + 0.93x(2). V.o: 9.467925e-03.