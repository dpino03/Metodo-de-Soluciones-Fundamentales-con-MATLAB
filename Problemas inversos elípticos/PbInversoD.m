function [a_opt] = PbInversoD(R,r,xc,dist,Nb,lambda0,lambda1,A,h,g,pert,a0) 
% -------------------------------------------------------------------------
% Consideramos un problema (P), de la forma:
%
% -Lap(u) + a*u = h en B(0;R)
% u = g             sobre fr(B(0;R)) 
%
% Sea u_a el valor de la solución del problema asociado al valor a = a(x).
%
% Consideraremos a = A(x) como el valor de a buscado. Dada la solución u_A
% asociada a la solución de (P_a) para este valor de a, realizaremos una 
% medición en D, la bola de centro xc y radio r, del valor de u_A.
%
% El espacio donde consideraremos a(x) serán los polinomios de grado 1.
%
% Entradas:
%
% R = Radio de la bola.
% r = Radio del dominio donde se realiza la medición.
% xc = Centro de la bola (vector 2x1).
% dist = Distancia entre los puntos del mallado regular interior para el
% cálculo de la solución particular u_P.
% Nb = Número de puntos sobre la frontera para la condición Dirichlet
% lambda0 = Factor para ubicar las fuentes fuera de la bola.
% lambda1 = Constante que aparece en la definición de las funciones
% radiales f y F, relacionadas por Lap(F) = f.
% A = Función que describe la EDP para la que realizaremos la medición.
% h = Función que describe la EDP.
% g = Función que describe las condiciones Dirichlet en la frontera.
% pert = Porcentaje de ruido en la medición. Por ejemplo, pert = 1.e-2 es
% una perturbación aleatoria del 1%.
% a0 = Aproximación de los coeficientes del polinomio buscado, para poder
% aplicar de forma eficiente fmincon.
%
% Salida:
%
% Valor del polinomio a(x) óptimo.
%
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Vamos a resolver en primer lugar el problema (P) para el valor a = A(x)
% sobre el que realizaremos la medición.
% -------------------------------------------------------------------------

[u_A] = MSF_MSPEliptDir(R,xc,dist,Nb,lambda0,lambda1,A,h,g);

% Definimos ahora la función eta(x) como el valor de u_A restringido a la
% bola B(xc;r). Como u_A está definido en todo B(xc;R), basta considerar
% eta(x) como u_A, y evaluarla en puntos de B(xc;r):

eta = @(x) u_A(x);

% Esta será nuestra medición en D.

% Vamos a considerar una perturbación de la medición, es decir,
% consideraremos la función eta modificada aleatoriamente, lo que 
% representa el ruido en la medición. Para considerar una medición exacta,
% basta considerar pert = 0. El parámetro pert determinará el porcentaje de
% ruido aleatorio.

ruido = (2*rand(size(eta([0,0])))-1);  % Ruido aleatorio
eta_pert = @(x) eta(x)+pert*ruido.*eta(x);

% Nótese que el ruido debe de ser fijo, y no depender de x, pues en caso 
% contrario la función eta_pert no cumpliría con la regularidad suficiente 
% para poder aplicar fmincon.

% -------------------------------------------------------------------------
% Como a(x) toma valores en los polinomios de grado menor o igual que uno,
% haremos una identificación del polinomio  
% a(x) = a_1 + a_2*x(1) + a_3*x(2) 
% con el vector coef_a = [a_1,a_2,a_3].
% -------------------------------------------------------------------------

% Función para evaluar el polinomio a(x):

function pol = polinomio(a,x)
    pol = a(1)+a(2)*x(1)+a(3)*x(2);
end

% -------------------------------------------------------------------------
% Definimos ahora el funcional a minimizar, J(a).
% -------------------------------------------------------------------------

function J_a = funcJ(a_coef)
    a_func = @(x) polinomio(a_coef,x);
    [u_a] = MSF_MSPEliptDir(R,xc,dist,Nb,lambda0,lambda1,a_func,h,g);
    f = @(x) norm(eta_pert(x)-u_a(x)).^2;
    % Para integrar, definimos el integrando con un cambio de variable a 
    % coordenadas polares.
    integ = @(rho,theta) arrayfun(@(r,t) f([xc(1)+r*cos(t);xc(2)+r*sin(t)]),rho,theta).*rho;
    J_a = 0.5*integral2(integ,0,r,0,2*pi);
end

% -------------------------------------------------------------------------
% Resolución con fmincon.
% -------------------------------------------------------------------------

% Llamamos a fmincon:

options = optimoptions('fmincon', ...
    'Algorithm', 'sqp', ...             % Algoritmo SQP (Sequential Quadratic Programming)
    'Display', 'iter', ...              % Información por pantalla en cada iteración 
    'OptimalityTolerance', 1e-5, ...    % Tolerancia para considerar que se ha alcanzado la optimalidad 
    'StepTolerance', 1e-6, ...          % Tolerancia para considerar que los pasos ya no cambian mucho
    'FiniteDifferenceType', 'central'); % Usa diferencias finitas centradas para estimar derivadas 

a_opt = fmincon(@funcJ,a0,[],[],[],[],[],[],[],options);

end





