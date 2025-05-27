function [a_opt] = PbInversogamma(R,xc,dist,Nb,lambda0,lambda1,A,h,g,pert,a0) 
% -------------------------------------------------------------------------
% Consideramos un problema (P), de la forma:
%
% -Lap(u) + a*u = h en B(0;R)
% u = g             sobre fr(B(0;R)) 
%
% Sea u_a el valor de la solución del problema asociado al valor a = a(x).
%
% Consideraremos a = A(x) como el valor de a buscado. Dada la solución u_A
% asociada a la solución de (P) para este valor de a, realizaremos una 
% medición en gamma, la semicircunferencia superior de la frontera de
% B(xc;R), de la derivada normal de u_A.
%
% El espacio donde consideraremos a(x) serán los polinomios de grado 1.
%
% Entradas:
%
% R = Radio de la bola.
% xc = Centro de la bola (vector 2x1).
% dist = Distancia entre los puntos del mallado regular interior para el
% cálculo de la solución particular u_P.
% Nb = Número de puntos sobre la frontera para la condición Dirichlet
% lambda0 = Factor para ubicar las fuentes fuera de la bola.
% lambda1 = Constante que aparece en la definición de las funciones
% radiales f y F, relacionadas por ΔF = f.
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

% Definimos ahora la función eta(x) como el valor de de la derivada normal
% de u_A en la semicircunferencia superior de la frontera.

% Calculamos el gradiente de forma numérica.

h_df = 1e-5;  % Tamaño del paso para la diferencia finita

grad_u_A = @(x) [ ...
    (u_A(x+[h_df; 0])-u_A(x-[h_df; 0]))/(2*h_df); ...
    (u_A(x+[0; h_df])-u_A(x-[0; h_df]))/(2*h_df) ...
];

% El vector normal en una bola de centro xc y radio R viene dado por la
% siguiente expresión:

n = @(x) (x-xc)./R;

% Luego tenemos la siguiente aproximación numérica de la derivada normal:

eta = @(x) grad_u_A(x).*n(x);

% Vamos a considerar una perturbación de la medición, es decir,
% consideraremos la función eta modificada aleatoriamente, lo que 
% representa el ruido en la medición. Para considerar una medición exacta,
% basta considerar pert = 0. El parámetro pert determinará el porcentaje de
% ruido aleatorio.

ruido = (2*rand(size(eta([0,0])))-1);  % Ruido aleatorio
eta_pert = @(x) eta(x)+pert*ruido.*eta(x);


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
    grad_u_a = @(x) [ ...
    (u_a(x+[h_df; 0])-u_a(x - [h_df; 0]))/(2*h_df); ...
    (u_a(x+[0; h_df])-u_a(x - [0; h_df]))/(2*h_df) ...
    ];
    eta_a = @(x) grad_u_a(x).*n(x); 
    f = @(x) norm(eta_pert(x)-eta_a(x))^2;
    % Cambio a coordenadas polares para la semicircunferencia superior.
    integ = @(theta) arrayfun(@(t) f([xc(1)+R*cos(t);xc(2)+R*sin(t)])*R,theta);
    % Integral sobre theta en [0,pi] (semicircunferencia superior).
    J_a = 0.5*integral(integ,0,pi,'ArrayValued',true);
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