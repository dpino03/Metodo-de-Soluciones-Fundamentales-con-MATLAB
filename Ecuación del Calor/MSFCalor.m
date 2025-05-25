function [u] = MSFCalor(R,xc,T,M,u0,dist,Nb,lambda0,lambda1,h,g,u_exacta)
% -------------------------------------------------------------------------
% Consideramos un problema de la forma:
%
% u_t - Lap(u) = h(x,t) en B(xc;R)x(0,T]
% u = g(x,t)            sobre fr(B(xc;R))x(0,T]
% u(x,0) = u0           en B(xc;R)x{0}
%
% Entradas:
%
% R = Radio de la bola.
% xc = Centro de la bola (vector 2x1).
% T = Tiempo final.
% M = Número de subintervalos en los que dividiremos [0,T].
% u0 = Valor de la función en el instante inicial t = 0.
% dist = Distancia entre los puntos del mallado regular interior para el
% cálculo de la solución particular u_P en cada problema elíptico. Se usará
% al invocar la función MSF_MSPEliptDir.
% Nb = Número de puntos de la frontera para aplicar el MSF en cada problema
% elíptico. Se usará al invocar la función MSF_MSPEliptDir.
% lambda0 = Factor para ubicar las fuentes fuera de la bola. Se usará al 
% invocar la función MSF_MSPEliptDir.
% lambda1 = Constante que aparece en la definición de las funciones
% radiales f y F, relacionadas por ΔF = f. Se usará al invocar la función 
% MSF_MSPEliptDir.
% h = Función que describe la EDP.
% g = Función que describe las condiciones de Dirichlet en la frontera para
% cada tiempo t.
% u_exacta = Solución exacta para realizar la comparación.
%
% Salida:
%
% Aproximación de la solución u(x) por el MSF-MSP y representación gráfica.
% Comparación de la solución del método con la solución exacta en cada
% tiempo.
%
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% En primer lugar, vamos a realizar la partición uniforme del intervalo
% [0,T] en los N subintervalos, es decir, tomar N+1 puntos.
% -------------------------------------------------------------------------

t = linspace(0,T,M+1); % Dividimos [0,T] en N+1 puntos
k = T/M;               % Este es el paso, pues la partición es uniforme

% -------------------------------------------------------------------------
% Definimos ahora un bucle que nos va a devolver la solución en cada
% instante de tiempo, resolviendo el problema elíptico visto en la parte
% teórica en cada instante de tiempo.
% -------------------------------------------------------------------------

% Definimos u como una celda que irá almacenando en cada entrada  la 
% solución en el instante de tiempo correspondiente como función anónima.

u = cell(1,M+1); 

% Inicializamos con la solución en el instante inicial t=0.

u{1} = u0;    

% Definimos el siguiente valor como escalar y como función anónima, ya que
% en diferentes contextos en el bucle lo necesitaremos de una forma u otra.

k1 = 1/k; 
a = @(x) k1;

% --- Desarrollo del bucle ---

for n = 1:M
    un = u{n};  
    tn1 = t(n+1); 
    h1 = @(x) h(x,tn1);
    h_elipt = @(x) h1(x) + k1*un(x);
    g1 = @(x) g(x,tn1);
    [uu] = MSF_MSPEliptDir(R,xc,dist,Nb,lambda0,lambda1,a,h_elipt,g1);
    u{n+1} = uu;
end

% -------------------------------------------------------------------------
% Seguiremos con la representación gráfica de la solución. Lo haremos de
% forma que veamos cómo varía la solución a lo largo del tiempo.
% -------------------------------------------------------------------------

% --- Construcción del mallado ---

% Número de puntos de evaluación en el radio y en el ángulo.

Nr = 70;
Ntheta = 70;

% Creamos una malla en coordenadas polares dentro de la bola.

r = linspace(0,R,Nr);                 % Radios desde 0 hasta R
theta_mesh = linspace(0,2*pi,Ntheta); % Ángulos de 0 a 2pi
[Rad,Theta] = meshgrid(r,theta_mesh); % Malla en coordenadas polares

% Convertimos la malla a coordenadas cartesianas.

X = Rad.*cos(Theta) + xc(1);
Y = Rad.*sin(Theta) + xc(2);

% -------------------------------------------------------------------------
% Vamos a ir construyendo a la vez la representación de la solución, de la
% solución exacta y del error, ya que el proceso es análogo. Es por esto, 
% para no doblar los bucles, haremos los desarrollos simultáneamente.
% -------------------------------------------------------------------------

% Inicializamos la solución, el error y la solución exacta.
  
U = nan(size(X)); 
E = nan(size(X)); 
Ue = nan(size(X));

% Graficamos la solución, solución exacta y el error.

for n = 1:M+1
    for i = 1:size(X,1)
        for j = 1:size(X,2)
            U(i,j) = u{n}([X(i,j);Y(i,j)]);                      % Solución en cada punto de la malla
            E(i,j) = abs(U(i,j)-u_exacta([X(i,j);Y(i,j)],t(n))); % Error en cada punto de la malla
            Ue(i,j) = u_exacta([X(i,j);Y(i,j)],t(n));            % Solución exacta en cada punto de la malla
        end
    end

    % Gráfica de la la solución y el error para el tiempo t(n).

    figure(1)
    surf(X,Y,U,'EdgeColor','none');                       % Gráfica 3D sin bordes
    colorbar;                                             % Barra de color
    title(['Solución del método en t = ',num2str(t(n))]); % Título
    xlabel('x_1');ylabel('x_2');zlabel('u(x)');           % Etiquetas de los ejes
    axis('equal','tight');                                % Configuración de los ejes
    view(3);                                              % Vista en 3D
    shading('interp');                                    % Sombreado por interpolación
    pause(0.1)                                            % Pausa de 0.1 segundos 
   
    figure(2)
    surf(X,Y,Ue,'EdgeColor','none'); 
    colorbar; 
    title(['Solución exacta en t = ',num2str(t(n))]); 
    xlabel('x_1');ylabel('x_2');zlabel('u(x)'); 
    axis('equal','tight'); 
    view(3); 
    shading('interp'); 
    pause(0.1); 
    
    figure(3)
    surf(X,Y,E,'EdgeColor','none'); 
    colorbar; 
    title(['Error absoluto en t = ',num2str(t(n))]); 
    xlabel('x_1');ylabel('x_2');zlabel('u(x)'); 
    axis('equal','tight'); 
    view(3); 
    shading('interp'); 
    pause(0.1)
end

saveas(figure(1),'SolucionMetPbCalorDir.png');
saveas(figure(2),'SolucionExactaPbCalorDir.png');
saveas(figure(3),'ErrorPbCalorDir.png');

end