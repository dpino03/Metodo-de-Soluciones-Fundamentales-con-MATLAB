function [u] = MSFCalorDir(R,xc,T,M,u0,dist,Nb,lambda0,lambda1,h,g)
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
%
% Salida:
%
% Aproximación de la solución u(x) por el MSF-MSP.
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

end