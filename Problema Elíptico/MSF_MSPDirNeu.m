function [u] = MSF_MSPDirNeu(R,xc,dist,Nb1,Nb2,lambda0,lambda1,a,h,g,xi,u_exacta)
% -------------------------------------------------------------------------
% Consideramos un problema de la forma:
%
% -Lap(u) + a*u = h en B(0;R)
% u = g             sobre Gamma1 contenida en la frontera de B(0;R)
% grad(u)*n = xi    sobre Gamma2 contenida en la frontera de B(0;R)
%
% Tomaremos como Gamma1 la semicircunferencia superior de la frontera de la
% bola y como Gamma2 la semicircunferencia inferior.
%
% Entradas:
%
% R = Radio de la bola.
% xc = Centro de la bola (vector 2x1).
% dist = Distancia entre los puntos del mallado regular interior para el
% cálculo de la solución particular u_P.
% Nb1 = Número de puntos Gamma1 para la condición Dirichlet
% Nbgamma = Número de puntos sobre Gamma2 para imponer la condición de tipo
% Neumann.
% lambda0 = Factor para ubicar las fuentes fuera de la bola.
% lambda1 = Constante que aparece en la definición de las funciones
% radiales f y F, relacionadas por ΔF = f.
% a = Función que describe la EDP.
% h = Función que describe la EDP.
% g = Función que describe las condiciones Dirichlet en la frontera.
% xi = Función que describe las condiciones Neumann en la frontera.
% u_exacta = Solución exacta para realizar la comparación
%
% Salida:
%
% Aproximación de la solución u(x) por el MFS-MSP y representación gráfica.
% Comparación de la solución del método con la solución exacta.
%
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
%  En primer lugar, tomaremos una partición de 0 a pi con Nb1 puntos para 
%  tener Nb1 puntos equidistantes sobre Gamma1 en la frontera de la bola.
% -------------------------------------------------------------------------

% ----- Puntos sobre Gamma1 (semicircunferencia superior) -----

theta = linspace(0,pi,Nb1+1); % Dividimos la circunferencia en Nb1+1 puntos
theta(end) = [];              % Eliminamos el último para evitar el punto repetido

% Coordenadas de los puntos frontera.

m1 = R*cos(theta)+xc(1); % Coordenadas x1 de la frontera
m2 = R*sin(theta)+xc(2); % Coordenadas x2 de la frontera
puntos_Gamma1 = [m1;m2]; % Matriz 2xNb1 de puntos sobre la frontera

% -------------------------------------------------------------------------
% A continuación, elegimos los puntos fuente fuera de la bola para Gamma1,
% en una circunferencia de centro xc y radio R+lambda0.
% -------------------------------------------------------------------------

% ----- Puntos fuente fuera de la bola para Gamma1 -----

% Las fuentes están en una circunferencia de radio R*lambda0.

s1_G1 = lambda0*R*cos(theta)+xc(1); % Coordenadas x1 de las fuentes
s2_G1 = lambda0*R*sin(theta)+xc(2); % Coordenadas x2 de las fuentes
fuentes_Gamma1 = [s1_G1;s2_G1];     % Matriz 2xNb1 de puntos fuentes fuera de la bola

% -------------------------------------------------------------------------
% Vamos a elegir ahora los Nb2 puntos sobre Gamma2, la
% semicircunferencia inferior, para luego imponer la condición de tipo
% Neumann.
% -------------------------------------------------------------------------

% ----- Puntos sobre Gamma2 (semicircunferencia inferior) -----

theta = linspace(pi,2*pi,Nb2+1); % Dividimos la semicircunferencia en Nb2+1 puntos
theta(end) = [];                 % Eliminamos el último para evitar el punto repetido

% Coordenadas de los puntos sobre Gamma2

z1 = R*cos(theta)+xc(1); % Coordenadas z1 de Gamma2
z2 = R*sin(theta)+xc(2); % Coordenadas z2 de Gamma2
puntos_Gamma2 = [z1;z2]; % Matriz 2xNb2 de puntos sobre Gamma2

% -------------------------------------------------------------------------
% Elegimos los puntos fuente fuera de la bola para Gamma2, en una
% circunferencia de centro xc y radio R+lambda0.
% -------------------------------------------------------------------------

% ----- Puntos fuente fuera de la bola para Gamma1 -----

% Las fuentes están en una circunferencia de radio R*lambda0.

s1_G2 = lambda0*R*cos(theta)+xc(1); % Coordenadas x1 de las fuentes
s2_G2 = lambda0*R*sin(theta)+xc(2); % Coordenadas x2 de las fuentes
fuentes_Gamma2 = [s1_G2;s2_G2];     % Matriz 2xNb1 de puntos fuentes fuera de la bola

% -------------------------------------------------------------------------
% Definimos el total de las fuentes, las de Gamma1 y las de Gamma2.
% -------------------------------------------------------------------------

puntos_fuentes = [fuentes_Gamma1,fuentes_Gamma2]; % Matriz 2xNb

% -------------------------------------------------------------------------
% Finalmente, consideramos el mallado interior con los puntos de campo.
% -------------------------------------------------------------------------

% Creamos una malla rectangular que cubra el área de la bola.

x = (xc(1)-R):dist:(xc(1)+R);
y = (xc(2)-R):dist:(xc(2)+R);
[X,Y] = meshgrid(x,y);

% Filtramos los puntos estrictamente dentro del círculo de radio R.

distancias = sqrt((X-xc(1)).^2 + (Y-xc(2)).^2);
dentro_circ = distancias < R;                   % Condición estricta

% Coordenadas de los puntos interiores

puntos_campo_x = X(dentro_circ);
puntos_campo_x = puntos_campo_x'; % Como vector fila
puntos_campo_y = Y(dentro_circ);
puntos_campo_y = puntos_campo_y'; % Como vector fila

puntos_campo = [puntos_campo_x;puntos_campo_y];

Nf = length(puntos_campo_x);

% -------------------------------------------------------------------------
% Ahora, definiremos la solución fundamental del Laplaciano, así como las 
% funciones f y F.
% -------------------------------------------------------------------------

% Añadimos valores de tolerancia en las definiciones para evitar
% singularidades.

G = @(r) -1/(2*pi)*log(max(r,1e-10));

f = @(r) (r <= lambda1).*(1 - max(r,1e-10)/lambda1).^2;

F = @(r) (r <= lambda1) .* (r^4/(16*lambda1^2) - (2*r^3)/(9*lambda1) + r^2/4) + ...
         (r > lambda1) .* (13*lambda1^2 + lambda1^2/12*log(max(r/lambda1,1e-10)));

% Vamos a definir también las derivadas de F y G como funciones de r, que
% es lo que necesitaremos para la derivada normal de estas funciones como
% se ha visto en el desarrollo teórico.

Fprim = @(r) (r<=lambda1) * (r^3/(4*lambda1^2)-(2*r^2)/(3*lambda1)+r/2) + ...
             (r>lambda1) * (lambda1^2/(12*max(r,1e-10)));

Gprim = @(r) -1./(2*pi*max(r,1e-10));

% -------------------------------------------------------------------------
% Procedemos con la definición de la matriz A que define el sistema a
% resolver.
% -------------------------------------------------------------------------

% ----- Construcción de la matriz A -----

% Vamos a construir la matriz por bloques:

Nb = Nb1+Nb2;

C = zeros(Nf,Nf);

for i = 1:Nf
    r_i = puntos_campo(:,i);
    for j = 1:Nf
        r_j = puntos_campo(:,j);
        C(i,j) = -f(norm(r_i-r_j)) + a(r_i)*F(norm(r_i-r_j));
    end
end

D = zeros(Nf,Nb);

for i = 1:Nf
    r_i = puntos_campo(:,i);
    for j = 1:Nb
        s_j = puntos_fuentes(:,j);
        D(i,j) = a(r_i)*G(norm(r_i-s_j));
    end
end

E = zeros(Nb1,Nf);

for i = 1:Nb1
    x_i = puntos_Gamma1(:,i);
    for j = 1:Nf
        r_j = puntos_campo(:,j);
        E(i,j) = F(norm(x_i-r_j));
    end
end

L = zeros(Nb1,Nb);

for i = 1:Nb1
    x_i = puntos_Gamma1(:,i);
    for j = 1:Nb
        s_j = puntos_fuentes(:,j);
        L(i,j) = G(norm(x_i-s_j));
    end
end

P = zeros(Nb2,Nf);

for i = 1:Nb2
    z_i = puntos_Gamma2(:,i);
    for j = 1:Nf
        r_j = puntos_campo(:,j);
        P(i,j) = Fprim(norm(z_i-r_j));
    end
end

Q = zeros(Nb2,Nb);

for i = 1:Nb2
    z_i = puntos_Gamma2(:,i);
    for j = 1:Nb
        s_j = puntos_fuentes(:,j);
        Q(i,j) = Gprim(norm(z_i-s_j));
    end
end

% Finalmente obtenemos la matriz A:

A = [C,D;E,L;P,Q]; % Matriz (Nf+Nb)x(Nf+Nb)

% La matriz A está mal condicionada para valores de N grandes, podemos
% realizar distintos métodos de precondicionamiento para subsanar este
% problema.

% Consideraremos la siguiente regularización:

lambda = 1e-6;               % Valor pequeño para regularización
A = A + lambda*eye(size(A)); 

% -------------------------------------------------------------------------
%  Creamos un vector con los términos independientes del sistema a resolver
% -------------------------------------------------------------------------

% Lo haremos igualmente por bloques, uno correspondiente a la EDP en los
% puntos de campo, otro a la condición de frontera de tipo Dirichlet sobre 
% Gamma1 y otro a la tipo Neumann sobre Gamma2.

b1 = zeros(1,Nf);

for i = 1:Nf
    r_i = puntos_campo(:,i);
    b1(i) = h(r_i);
end

b2 = zeros(1,Nb1);

for i = 1:Nb1
    m_i = puntos_Gamma1(:,i);
    b2(i) = g(m_i);
end

b3 = zeros(1,Nb2);

for i = 1:Nb2
    z_i = puntos_Gamma2(:,i);
    b3(i) = xi(z_i);
end

b = [b1,b2,b3];

b = b'; % Lo ponemos como vector columna de Nf+Nb componentes.

% -------------------------------------------------------------------------
% Resolvemos el sistema lineal dado por A*coefs = b, de modo que obtenemos
% los coeficientes alpha_j y beta_k necesarios para poder defininir la 
% solución. Construimos dicha solución.
% -------------------------------------------------------------------------

% ----- Resolver el sistema para los coeficientes -----

% Resolución como sistema lineal:

coefs = A\b;

% ----- Construcción de la función -----

% Una vez tenemos los coeficientes, podemos definir la función u = u(x):

beta = coefs(1:Nf);
alpha = coefs(Nf+1:end);

% Definimos u con una operación vectorial.

u = @(x) dot(beta,arrayfun(@(j) F(norm(x-puntos_campo(:,j))),1:size(puntos_campo,2))) + ...
         dot(alpha,arrayfun(@(k) G(norm(x-puntos_fuentes(:,k))),1:size(puntos_fuentes,2)));

% -------------------------------------------------------------------------
%  Vamos a representar la solución obtenida a través de un mallado de la
%  bola de centro xc y radio R. Para ello tomaremos una partición del radio
%  y otra del ángulo, para poder calcular así un mallado en coordenadas
%  polares. Una vez tengamos este mallado, lo pasamos a coordenadas
%  cartesianas de nuevo, obteniendo así las matrices X1 y X2.
%  Posteriormente, construimos la solución del método sobre este mallado. 
%  Como la función está definida por recursión, es más óptimo construirla 
%  paso apaso en el mallado concreto.
% -------------------------------------------------------------------------

% ----- Evaluación de la solución aproximada en un mallado dentro de la bola -----

Nr = 200;     % Número de puntos de evaluación en el radio
Ntheta = 200; % Número de puntos de evaluación en el ángulo

% Creamos una malla en coordenadas polares dentro de la bola,

r = linspace(0,R,Nr);                 % Radios desde 0 hasta R
theta_mesh = linspace(0,2*pi,Ntheta); % Ángulos de 0 a 2pi
[Rad,Theta] = meshgrid(r,theta_mesh); % Crear la malla en coordenadas polares

% Convertimos la malla a coordenadas cartesianas.

X1 = Rad.*cos(Theta) + xc(1);
X2 = Rad.*sin(Theta) + xc(2);

% Evaluamos la solución aproximada u en cada punto de la malla.

U = zeros(size(X1));

% Usamos la función ya definida

for i = 1:Nr
    for j = 1:Ntheta
        x_eval = [X1(i,j);X2(i,j)];  % Punto en el que evaluamos u
        U(i,j) = u(x_eval);  
    end
end

% -------------------------------------------------------------------------
%  Evaluamos la solución exacta en el mallado.
% -------------------------------------------------------------------------

ue =  zeros(size(X1));

for i = 1:Nr
    for j = 1:Ntheta
        x_eval = [X1(i,j);X2(i,j)];  % Punto en el que evaluamos u_exacta
        ue(i,j) = u_exacta(x_eval);  
    end
end

% -------------------------------------------------------------------------
%  Graficamos la solución del método y exacta, y tomamos el error absoluto
%  y relativo en cada punto. El error que tomaremos como referencia será la
%  norma infinito del error relativo.
% -------------------------------------------------------------------------

% ----- Gráficas -----

figure(1)
surf(X1,X2,U);
xlabel('x_1');
ylabel('x_2');
zlabel('u(x_1,x_2)');
title('Solución con MSF-MSP')
shading interp; % Para suavizar la visualización de la superficie
colorbar;       % Añadir barra de colores

figure(2)
surf(X1,X2,ue);
xlabel('x_1');
ylabel('x_2');
zlabel('u(x_1,x_2)');
title('Solución exacta')
shading interp;  
colorbar;  

figure(3)
surf(X1,X2,ue-U);
xlabel('x_1');
ylabel('x_2');
zlabel('Error absoluto');
title('Error absoluto del MSF-MSP respecto de la solución exacta')
shading interp; 
colorbar; 

epsilon = 1e-10;
max_ue = max(abs(ue),[],'all'); % Valor máximo de la solución exacta en todo el dominio
error_relativo = abs(ue-U)./max(abs(ue),epsilon*max_ue);

figure(4)
surf(X1, X2,error_relativo);
xlabel('x_1');
ylabel('x_2');
zlabel('Error relativo');
title('Error relativo del MSF-MSP respecto de la solución exacta')
shading interp;
colorbar;

err_rel_max = max(max(abs((ue-U)./max(ue,1.e-10))));
err_abs_max = max(max(abs(ue-U)));

% Como es una matriz, hay que aplicar dos veces máximo.

disp(['El error relativo máximo es ',num2str(err_rel_max)]);
disp(['El error absoluto máximo es ',num2str(err_abs_max)]);

% -------------------------------------------------------------------------
% Graficar puntos del mallado, Gamma1, Gamma2 y fuentes.
% -------------------------------------------------------------------------

figure(5);
hold on;

% 1. Puntos del mallado.
scatter(puntos_campo(1,:), puntos_campo(2,:), 10, 'b', 'filled');

% 2. Puntos de Gamma1.
scatter(puntos_Gamma1(1,:), puntos_Gamma1(2,:), 30, 'r', 'filled');

% 3. Puntos de Gamma2.
scatter(puntos_Gamma2(1,:), puntos_Gamma2(2,:), 30, 'g', 'filled');

% 4. Fuentes fuera de la bola.

scatter(fuentes_Gamma1(1,:), fuentes_Gamma1(2,:), 50, 'm', 'x', 'LineWidth', 1.5);
scatter(fuentes_Gamma2(1,:), fuentes_Gamma2(2,:), 50, 'm', 'x', 'LineWidth', 1.5);

% Dibujar la frontera de la bola

theta_circ = linspace(0, 2*pi, 100);
plot(R*cos(theta_circ) + xc(1), R*sin(theta_circ) + xc(2), 'k-', 'LineWidth', 1.5);

% Configuración de la gráfica

xlabel('x_1', 'Color', 'k', 'FontWeight', 'bold');
ylabel('x_2', 'Color', 'k', 'FontWeight', 'bold');
title('Distribución de puntos y fuentes', 'Color', 'k', 'FontWeight', 'bold');
legend('Puntos interiores r_i ', 'Puntos m_i sobre \Gamma_1 (Dirichlet)', ...
        'Puntos z_i sobre \Gamma_2 (Neumann)', 'Fuentes', 'Location', 'best');
grid on;
axis equal;
hold off;


end