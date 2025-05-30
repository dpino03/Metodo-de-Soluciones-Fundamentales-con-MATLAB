function [u] = MSF_MSPEliptDir(R,xc,dist,Nb,lambda0,lambda1,a,h,g)
% -------------------------------------------------------------------------
% Consideramos un problema de la forma:
%
% -Δu + a*u = h en B(0;R)
% u = g         sobre fr(B(0;R))
%
% Entradas:
%
% R = Radio de la bola.
% xc = Centro de la bola (vector 2x1).
% dist = Distancia entre los puntos del mallado regular interior para el
% cálculo de la solución particular u_P.
% Nb = Número de puntos de la frontera para aplicar el MSF
% lambda0 = Factor para ubicar las fuentes fuera de la bola.
% lambda1 = Constante que aparece en la definición de las funciones
% radiales f y F, relacionadas por ΔF = f.
% a = Función que describe la EDP.
% h = Función que describe la EDP.
% g = Función que describe las condiciones de Dirichlet en la frontera.
% u_exacta = Solución exacta para realizar la comparación
%
% Salida:
%
% Aproximación de la solución u(x) por el MSF-MSP.
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
%  En primer lugar, tomaremos una partición de 0 a 2*pi con N puntos para 
%  tener Nb puntos equidistantes sobre la frontera de la bola.
% -------------------------------------------------------------------------

% ----- Puntos sobre la frontera (circunferencia de radio R) -----

theta = linspace(0,2*pi,Nb+1); % Dividimos la circunferencia en Nb+1 puntos
theta(end) = [];               % Eliminamos el último para evitar el punto repetido

% Coordenadas de los puntos frontera.

x1_front = R*cos(theta) + xc(1);       % Coordenadas x1 de la frontera
x2_front = R*sin(theta) + xc(2);       % Coordenadas x2 de la frontera
puntos_frontera = [x1_front;x2_front]; % Matriz 2xN de puntos sobre la frontera

% -------------------------------------------------------------------------
% A continuación, elegimos los puntos fuente fuera de la bola, en una
% circunferencia de centro xc y radio R+lambda0.
% -------------------------------------------------------------------------

% ----- Puntos fuente fuera de la bola -----

% Las fuentes están en una circunferencia de radio R*lambda0.

s1 = lambda0*R*cos(theta) + xc(1); % Coordenadas x1 de las fuentes
s2 = lambda0*R*sin(theta) + xc(2); % Coordenadas x2 de las fuentes
puntos_fuentes = [s1;s2];          % Matriz 2xNb de puntos fuentes fuera de la bola

% -------------------------------------------------------------------------
% Finalmente, consideramos el mallado interior con los puntos de campo
% -------------------------------------------------------------------------

% Creamos una malla rectangular que cubra el área de la bola.

x = (xc(1) - R):dist:(xc(1) + R);
y = (xc(2) - R):dist:(xc(2) + R);
[X,Y] = meshgrid(x,y);

% Filtramos los puntos estrictamente dentro del círculo de radio R.

distancias = sqrt((X - xc(1)).^2 + (Y - xc(2)).^2);
dentro_circ = distancias < R;                       % Condición estricta

% Coordenadas de los puntos interiores

puntos_campo_x = X(dentro_circ);
puntos_campo_x = puntos_campo_x'; % Como vector fila
puntos_campo_y = Y(dentro_circ);
puntos_campo_y = puntos_campo_y'; % Como vector fila

puntos_campo = [puntos_campo_x;puntos_campo_y];

Nf = length(puntos_campo_x);

% -------------------------------------------------------------------------
% Ahora, definiremos la solución fundamental del Laplaciano, así como las 
% funciones f y F
% -------------------------------------------------------------------------

% Añadimos valores de tolerancia en las definiciones para evitar
% singularidades.

G = @(r) -1/(2*pi)*log(max(r,1e-10));

f = @(r) (r <= lambda1).*(1 - max(r,1e-10)/lambda1).^2 + (r > lambda1).*0;

F = @(r) (r <= lambda1) .* (r^4/(16*lambda1^2) - (2*r^3)/(9*lambda1) + r^2/4) + ...
         (r > lambda1) .* (13*lambda1^2 + lambda1^2/12*log(max(r/lambda1,1e-10)));

% -------------------------------------------------------------------------
% Procedemos con la definición de la matriz A que define el sistema a
% resolver.
% -------------------------------------------------------------------------

% ----- Construcción de la matriz A -----

% Vamos a construir la matriz por bloques:

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
        D(i,j) =  a(r_i)*G(norm(r_i-s_j));
    end
end

E = zeros(Nb,Nf);

for i = 1:Nb
    x_i = puntos_frontera(:,i);
    for j = 1:Nf
        r_j = puntos_campo(:,j);
        E(i,j) = F(norm(x_i-r_j));
    end
end

L = zeros(Nb,Nb);

for i = 1:Nb
    x_i = puntos_frontera(:,i);
    for j = 1:Nb
        s_j = puntos_fuentes(:,j);
        L(i,j) = G(norm(x_i-s_j));
    end
end

% Finalmente obtenemos la matriz A:

A = [C,D;E,L]; % Matriz (Nf+Nb)x(Nf+Nb)

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
% puntos de campo, y otro a la condición de frontera de tipo Dirichlet.

b1 = zeros(1,Nf);

for i = 1:Nf
    r_i = puntos_campo(:,i);
    b1(i) = h(r_i);
end

b2 = zeros(1,Nb);

for i = 1:Nb
    x_i = puntos_frontera(:,i);
    b2(i) = g(x_i);
end

b = [b1,b2];

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

end