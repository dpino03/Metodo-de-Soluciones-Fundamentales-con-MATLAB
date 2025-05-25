function [u] = MSFLaplaceRectangulo(xc,l1,l2,N,lambda,g,u_exacta)
% -------------------------------------------------------------------------
% Resuelve el problema de Laplace en un rectángulo:
%
% Lap(u) = 0 en Omega = [xc(1)-l1/2,xc(1)+l1/2]x[xc(2)-l2/2,xc(2)+l2/2]
% u = g      sobre fr(Omega)
%
% Entradas:
% 
% xc = Centro del rectángulo (vector 2x1).
% l1 = Longitud de los lados horizontales (superior e inferior).
% l2 = Longitud de los lados verticales (izquierdo y derecho).
% N = Número de puntos de frontera (deberá ser múltiplo de 4).
% lambda = Factor de distancia para ubicar fuentes.
% g = Función de condición de Dirichlet en la frontera.
% u_exacta = Solución exacta.
%
% Salida:
%
% Aproximación de la solución u(x) por el MFS, error y representación 
% gráfica.
%
% -------------------------------------------------------------------------

% En primer lugar, aseguramos que tenemos un número de puntos total
% múltiplo de 4, de forma que podamos asegurar una distribución uniforme de
% los puntos sobre la frontera.

if mod(N,4) ~= 0
    error('N debe ser múltiplo de 4 para distribución uniforme en el rectángulo.');
end

% -------------------------------------------------------------------------
% Comenzamos con la generación de puntos en la frontera y de las fuentes.
% -------------------------------------------------------------------------

% Definimos el número de puntos en cada lado (N/4 por lado).

N_lado = N/4;

% Definimos las coordenadas de los vértices del rectángulo.

x_min = xc(1) - l1/2; 
x_max = xc(1) + l1/2;
y_min = xc(2) - l2/2; 
y_max = xc(2) + l2/2;

% Vamos ahora a definir cada uno de los lados del rectángulo, de forma que
% luego podamos trabajar con ellos.

% Lado inferior (de izquierda a derecha):
x_inf = linspace(x_min,x_max,N_lado+1); x_inf(end) = [];
y_inf = y_min*ones(1,N_lado);

% Lado derecho (de abajo a arriba):
y_der = linspace(y_min,y_max,N_lado+1); y_der(end) = [];
x_der = x_max * ones(1,N_lado);

% Lado superior (de derecha a izquierda):
x_sup = linspace(x_max, x_min, N_lado+1); x_sup(end) = [];
y_sup = y_max * ones(1,N_lado);

% Lado izquierdo (de arriba a abajo):
y_izq = linspace(y_max, y_min, N_lado+1); y_izq(end) = [];
x_izq = x_min * ones(1,N_lado);

% Combinamos ahora todos los puntos frontera.

puntos_frontera = [...
    [x_inf;y_inf], ...
    [x_der;y_der], ...
    [x_sup;y_sup], ...
    [x_izq;y_izq]];

% Definimos los vectores normales unitarios exteriores (ordenados según los
% lados).

normales = [...
    [zeros(1,N_lado); -ones(1,N_lado)], ...   % Lado inferior (normal hacia abajo)
    [ones(1,N_lado); zeros(1,N_lado)], ...    % Lado derecho (normal hacia la derecha)
    [zeros(1,N_lado); ones(1,N_lado)], ...    % Lado superior (normal hacia arriba)
    [-ones(1,N_lado); zeros(1,N_lado)]];      % Lado izquierdo (normal hacia la izquierda)

% Ya podemos definir la ubicación de las fuentes (desplazadas según la
% normal), siguiendo el parámetro lamda.

puntos_fuentes = puntos_frontera + lambda*normales;

% -------------------------------------------------------------------------
% Definimos la solución fundamental y la construcción del sistema lineal.
% -------------------------------------------------------------------------

% Introducimos la solución fundamental 2D del Laplaciano.

G = @(r) -1/(2*pi)*log(max(norm(r),1e-10));

% Construimos ahora la matriz del sistema.

A = zeros(N,N);

for i = 1:N
    x_i = puntos_frontera(:,i);
    for j = 1:N
        s_j = puntos_fuentes(:,j);
        A(i,j) = G(x_i - s_j);
    end
end

% Consideramos ahora el vector de condiciones de frontera

b = zeros(N,1);

for i = 1:N
    b(i) = g(puntos_frontera(:,i));
end

% -------------------------------------------------------------------------
% Pasamos a la resolución del sistema, con regularización. Podemos así
% definir la solución proporcionada por el método.
% -------------------------------------------------------------------------

At = transpose(A);
epsilon = 1e-8; % Parámetro de regularización
B = At*A + epsilon*eye(size(A));
d = At*b;
coefs = B\d;

% Definición de la solución aproximada

u = @(x) dot(coefs,arrayfun(@(j) G(x-puntos_fuentes(:,j)),1:N));

% -------------------------------------------------------------------------
% Procedemos con la evaluación en el dominio rectangular mediante un
% mallado.
% -------------------------------------------------------------------------

% Creamos en primer lugar la malla para evaluación.

Nmalla = 100;
x1_vals = linspace(x_min, x_max,Nmalla);
x2_vals = linspace(y_min, y_max,Nmalla);
[X1,X2] = meshgrid(x1_vals,x2_vals);

% Evaluamos la solución en la malla.

U = zeros(size(X1));

for i = 1:size(X1,1)
    for j = 1:size(X1,2)
        x_eval = [X1(i,j); X2(i,j)];
        % Sólo evaluamos si está dentro del rectángulo
        if x_eval(1) >= x_min && x_eval(1) <= x_max && ...
           x_eval(2) >= y_min && x_eval(2) <= y_max
            U(i,j) = u(x_eval);
        else
            U(i,j) = NaN;
        end
    end
end

% -------------------------------------------------------------------------
% Finalmente, representamos gráficamente todo lo anterior.
% -------------------------------------------------------------------------

figure(1)
surf(X1, X2, U);
xlabel('x_1');
ylabel('x_2');
zlabel('u(x_1,x_2)');
title(['Solución con MSF del Laplaciano en rectángulo [', ...
       num2str(x_min),',',num2str(x_max),']x[', ...
       num2str(y_min),',',num2str(y_max),']'])
shading interp;
colorbar;

% Ahora, vamos a visualizar la solución exacta.

U_exact = zeros(size(X1));

for i = 1:size(X1,1)
    for j = 1:size(X1,2)
        x_eval = [X1(i,j); X2(i,j)];
        if x_eval(1) >= x_min && x_eval(1) <= x_max && ...
           x_eval(2) >= y_min && x_eval(2) <= y_max
            U_exact(i,j) = u_exacta(x_eval);
         else
            U_exact(i,j) = NaN;
         end
    end
end
    
figure(2)   
surf(X1,X2,U_exact);  
xlabel('x_1');   
ylabel('x_2');  
zlabel('u(x_1,x_2)');  
title('Solución exacta');   
shading interp;    
colorbar;
    
% Finalmente, calculamos y representamos el error absoluto.
   
figure(3)
surf(X1, X2,abs(U-U_exact)); 
xlabel('x_1'); 
ylabel('x_2');  
zlabel('Error absoluto');
title('Error absoluto del MSF');
shading interp;
colorbar;


end

