function [u] = MSFLaplaceBola(R,xc,N,lambda,g,u_exacta)
% -------------------------------------------------------------------------
% Consideramos un problema de la forma:
%
% Lap(u) = 0 en B(xc;R)
% u = g      sobre la frontera de B(xc;R)
%
% Entradas:
% 
% R = Radio de la bola.
% xc = Centro de la bola (vector 2x1).
% N = Número de puntos de frontera para aplicar las condiciones de contorno.
% lambda = Factor para ubicar las fuentes fuera de la bola.
% g = Función que describe las condiciones de Dirichlet en la frontera.
% u_exacta = Solución exacta del problema.
%
% Salida:
%
% Aproximación de la solución u(x) por el MFS, error y representación 
% gráfica.
%
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
%  En primer lugar, tomaremos una partición de 0 a 2*pi con N puntos para 
%  tener N puntos equidistantes sobre la frontera de la bola.
% -------------------------------------------------------------------------

% ----- Puntos sobre la frontera (circunferencia de radio R) -----

theta = linspace(0,2*pi,N+1); % Dividimos la circunferencia en N+1 puntos
theta(end) = [];              % Eliminamos el último para evitar el punto repetido

% Coordenadas de los puntos frontera

x1_front = R*cos(theta) + xc(1);       % Coordenadas x1 de la frontera
x2_front = R*sin(theta) + xc(2);       % Coordenadas x2 de la frontera
puntos_frontera = [x1_front;x2_front]; % Matriz 2xN de puntos sobre la frontera

% -------------------------------------------------------------------------
% A continuación, elegimos los puntos fuente fuera de la bola, en una
% circunferencia de centro xc y radio R+lambda.
% -------------------------------------------------------------------------

% ----- Puntos fuente fuera de la bola -----

% Las fuentes están en una circunferencia de radio R*lambda

s1 = lambda*R*cos(theta) + xc(1); % Coordenadas x1 de las fuentes
s2 = lambda*R*sin(theta) + xc(2); % Coordenadas x2 de las fuentes
puntos_fuentes = [s1;s2];         % Matriz 2xN de fuentes fuera de la bola

% -------------------------------------------------------------------------
%  Ahora, definiremos la solución fundamental del Laplaciano, y
%  construiremos la matriz A = (a_{ij}), donde a_{i,j} = phi(x_i,s_j), para
%  x_i el i-ésimo punto tomado sobre la frontera y s_j el j-ésimo punto
%  fuente tomado sobre la frontera virtual.
% -------------------------------------------------------------------------

% ----- Definición de la solución fundamental -----

G = @(r) -1/(2*pi)*log(max(norm(r),1e-10)); % Solución fundamental 2D del Laplaciano

% Tomamos máximo en el argumento del logaritmo para evitar errores debidos
% a la singularidad.

% ----- Construcción de la matriz A -----

A = zeros(N,N);  % Matriz A de tamaño NxN

for i = 1:N
    x_i = puntos_frontera(:,i);    % Punto en la frontera x_i
    for j = 1:N
        s_j = puntos_fuentes(:,j); % Fuente s_j
        A(i,j) = G(x_i-s_j);       % Evaluación de la solución fundamental
    end
end

% -------------------------------------------------------------------------
%  Creamos un vector con la condición de frontera evaluada en cada punto x_i
%  que hemos tomado en la frontera
% -------------------------------------------------------------------------

% ----- Condiciones en la frontera -----

% Evaluamos la función g en los puntos de la frontera

b = zeros(N,1);

for i = 1:N
    b(i) = g(puntos_frontera(:,i)); % Condiciones de Dirichlet
end

% -------------------------------------------------------------------------
%  Resolvemos el sistema lineal dado por A*c = G, de modo que obtenemos los
%  coeficientes c_j necesarios para poder defininir la solución del MSF.
%  Construimos dicha solución.
% -------------------------------------------------------------------------

% ----- Resolver el sistema para los coeficientes -----

% La matriz A está mal condicionada para valores de N grandes, podemos
% realizar distintos métodos de precondicionamiento para subsanar este
% problema. Procederemos con una regularización de Tikhonov:

% Pasaremos del sistema Ac=b a Bc=d, donde B y b se construyen como sigue.

At = transpose(A);
epsilon = 1e-8;                  % Valor pequeño para regularización
B = At*A + epsilon*eye(size(A));
d = At*b;

coefs = B\d;                     % Resolver el sistema 


% ----- Función dada por el MSF -----

% Definimos la solución del método con una operación vectorial

u = @(x) dot(coefs,arrayfun(@(j) G(x-puntos_fuentes(:,j)),1:N));

% La función u_N viene dada por el producto escalar del vector de
% coeficientes por un vector de componentes phi(x,s_j).

% -------------------------------------------------------------------------
%  Vamos a representar la solución obtenida a través de un mallado de la
%  bola de centro xc y radio R. Para ello tomaremos una partición del radio
%  y otra del ángulo, para poder calcular así un mallado en coordenadas
%  polares. Una vez tengamos este mallado, lo pasamos a coordenadas
%  cartesianas de nuevo, obteniendo así las matrices X1 y X2.
%  Posteriormente, construimos la solución del MSF sobre este mallado. Como
%  la función está definida por recursión, es más óptimo construirla paso a
%  paso en el mallado concreto.
% -------------------------------------------------------------------------

% ----- Evaluación de u en un mallado dentro de la bola -----

Nr = 200;     % Número de puntos de evaluación en el radio
Ntheta = 200; % Número de puntos de evaluación en el ángulo

% Crear una malla en coordenadas polares dentro de la bola

r = linspace(0,R,Nr);                 % Radios desde 0 hasta R
theta_mesh = linspace(0,2*pi,Ntheta); % Ángulos de 0 a 2pi
[Rad,Theta] = meshgrid(r,theta_mesh); % Malla en coordenadas polares

% Convertir la malla a coordenadas cartesianas

X1 = Rad.*cos(Theta) + xc(1);
X2 = Rad.*sin(Theta) + xc(2);

% Evaluamos la solución aproximada u en cada punto de la malla

U = zeros(Nr,Ntheta);
for i = 1:Nr
    for j = 1:Ntheta
        x_eval = [X1(i,j); X2(i,j)]; % Punto en el que evaluamos u
        U(i,j) = u(x_eval);
    end
end

% -------------------------------------------------------------------------
%  Finalmente, evaluamos la solución exacta dada en los puntos de la malla,
%  de forma que podamos realizar la comparación.
% -------------------------------------------------------------------------

ue = zeros(Nr,Ntheta);
for i = 1:Nr
    for j = 1:Ntheta
        x_eval = [X1(i,j); X2(i,j)]; % Punto en el que evaluamos u_exacta
        ue(i,j) = u_exacta(x_eval);
    end
end

% -------------------------------------------------------------------------
%  Graficamos la solución del MSF y exacta, y tomamos el error absoluto y 
%  relativo en cada punto.
% -------------------------------------------------------------------------

% ----- Gráficas -----

figure(1)
surf(X1,X2,U);
xlabel('x_1');
ylabel('x_2');
zlabel('u(x_1,x_2)');
title(['Solución con MSF del Laplaciano con condiciones de tipo Dirichlet en B([', ...
        num2str(xc(1)),';',num2str(xc(2)),'],',num2str(R),')'])
shading interp; % Para suavizar la visualización de la superficie
colorbar;       % Añadir barra de colores

figure(2)
surf(X1,X2,ue);
xlabel('x_1');
ylabel('x_2');
zlabel('u(x_1,x_2)');
title(['Solución exacta Laplaciano con condiciones de tipo Dirichlet en B([', ...
        num2str(xc(1)),';',num2str(xc(2)),'],',num2str(R),')'])
shading interp;  
colorbar;  

figure(3)
surf(X1,X2,ue-U);
xlabel('x_1');
ylabel('x_2');
zlabel('Error absoluto');
title('Error absoluto del MSF respecto de la solución exacta')
shading interp; 
colorbar; 

epsilon = 1e-10;
max_ue = max(abs(ue),[],'all');  % Valor máximo de la solución exacta en todo el dominio
error_relativo = abs(ue-U)./max(abs(ue),epsilon*max_ue);

figure(4)
surf(X1,X2,error_relativo);
xlabel('x_1');
ylabel('x_2');
zlabel('Error relativo');
title('Error relativo del MSF respecto de la solución exacta')
shading interp;
colorbar;

err_rel_max = max(max(abs((ue-U)./max(ue,1.e-10))));
err_abs_max = max(max(abs(ue-U)));

% Como es una matriz, hay que aplicar dos veces máximo.

disp(['El error relativo máximo es ',num2str(err_rel_max),'.']);
disp(['El error absoluto máximo es ',num2str(err_abs_max),'.']);

end