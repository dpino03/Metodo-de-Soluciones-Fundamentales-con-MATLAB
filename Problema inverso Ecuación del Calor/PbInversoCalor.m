function [a_opt] = PbInversoCalor(R,xc,T,M,u0,dist,Nb,lambda0,lambda1,H,g,pert,h0) 
% -------------------------------------------------------------------------
% Consideramos un problema (P), de la forma:
%
% u_t - Lap(u) = h(x,t) en B(0;R)
% u = g                 sobre fr(B(0;R))x(0,T)
% u(x,0) = u0           en B(0;R)
%
% Sea u_h el valor de la solución del problema asociado al valor h = h(x).
%
% Consideraremos h = H(x) como el valor de a buscado. Dada la solución u_H
% asociada a la solución de (P) para este valor de h, realizaremos una 
% medición en gamma, la semicircunferencia superior de la frontera de
% B(xc;R), de la derivada normal de u_H en cada tiempo.
%
% El espacio donde consideraremos h(x) serán las constantes.
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
% H = Función que describe la EDP para la que realizaremos la medición.
% g = Función que describe las condiciones Dirichlet en la frontera.
% pert = Porcentaje de ruido en la medición. Por ejemplo, pert = 1.e-2 es
% una perturbación aleatoria del 1%.
% h0 = Aproximación de los coeficientes del polinomio buscado, para poder
% aplicar de forma eficiente fmincon.
%
% Salida:
%
% Valor del polinomio h(x) óptimo.
%
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Vamos a resolver en primer lugar el problema (P) para el valor h = AH(x)
% sobre el que realizaremos la medición.
% -------------------------------------------------------------------------

[u_H] = MSFCalorDir(R,xc,T,M,u0,dist,Nb,lambda0,lambda1,H,g);

% Definimos ahora la función nu(x,t) como el valor de de la derivada normal
% de u_H en la semicircunferencia superior de la frontera, para cada
% instante de tiempo

% Calculamos el gradiente y derivada normal para cada t_n

h_df = 1e-5;
grad_u_H = cell(1,M+1);
eta = cell(1,M+1);
eta_pert = cell(1,M+1);
normal = @(x) (x-xc)/R;    % Normal unitaria

% Hacemos un bucle en tiempo con el que calcularemos el gradiente, la
% derivada normal, y la perturbación de la medición.

for n = 1:M+1
    % Calculamos numéricamente el gradiente en t_n.

    grad_u_H{n} = @(x) [
        (u_H{n}(x+[h_df;0])-u_H{n}(x-[h_df;0]))/(2*h_df);
        (u_H{n}(x+[0;h_df])-u_H{n}(x-[0;h_df]))/(2*h_df);
    ];

    % Calculamos la derivada normal en t_n

    eta{n} = @(x) grad_u_H{n}(x).*normal(x);

    % Vamos a considerar una perturbación de la medición, es decir,
    % consideraremos la función nu modificada aleatoriamente en cada tiempo,lo 
    % que representa el ruido en la medición. Para considerar una medición 
    % exacta,basta considerar pert = 0. El parámetro pert determniará el 
    % porcentaje de ruido aleatorio.

    ruido = (2*rand(size(eta{n}([0,0])))-1);             % Ruido aleatorio
    eta_pert{n} = @(x) eta{n}(x)+pert*ruido.*eta{n}(x);
end

% -------------------------------------------------------------------------
% Definimos ahora el funcional a minimizar, J(h).
% -------------------------------------------------------------------------

function J_h = funcJ(h_coef)
    h_func = @(x,t) h_coef;
    [u_h] = MSFCalorDir(R,xc,T,M,u0,dist,Nb,lambda0,lambda1,h_func,g);
    % Inicializamos todo lo necesario.
    grad_u_h = cell(1,M+1);
    eta_h = cell(1,M+1);
    f = cell(1,M+1);
    integ = cell(1,M+1);
    J_h1 = zeros(1,M+1);
    % Hacemos el bucle en tiempo.
    for n = 1:M+1
        grad_u_h{n} = @(x) [
        (u_h{n}(x+[h_df;0])-u_h{n}(x-[h_df;0]))/(2*h_df);
        (u_h{n}(x+[0;h_df])-u_h{n}(x-[0;h_df]))/(2*h_df);
        ];
        % Definimos nu_h para cada tiempo t_n.
        eta_h{n} = @(x) grad_u_h{n}(x).*normal(x);
        % Definimos la función que integraremos para cada tiempo.
        f{n} = @(x) norm(eta_pert{n}(x)-eta_h{n}(x))^2;
        % Cambio a coordenadas polares para la semicircunferencia superior.
        integ{n} = @(theta) arrayfun(@(t) f{n}([xc(1)+R*cos(t);xc(2)+R*sin(t)])*R,theta);
        % Integral sobre theta en [0,pi] (semicircunferencia superior).
        J_h1(n) = 0.5*integral(integ{n},0,pi,'ArrayValued',true);
    end
    % Lo que minizaremos será la suma de las integrales en cada tiempo.
    J_h = sum(J_h1);
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

a_opt = fmincon(@funcJ,h0,[],[],[],[],[],[],[],options);

end