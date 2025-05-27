function [u] = MSFCalorGeneral(x1_eq,x2_eq,xc,T,M,u0,dist,Nb,lambda0,lambda1,h,g,u_exacta)
% -------------------------------------------------------------------------
% Consideramos un problema de la forma:
%
% u_t - Lap(u) = h(x,t) en B(xc;R)x(0,T]
% u = g(x,t)            sobre fr(B(xc;R))x(0,T]
% u(x,0) = u0           en B(xc;R)
%
% Entradas:
%
% x1_eq, x2_eq: ecuaciones paramétricas de la frontera.
% xc = Centro de la bola (vector 2x1).
% T = Tiempo final.
% M = Número de subintervalos en los que dividiremos [0,T].
% u0 = Valor de la función en el instante inicial t = 0.
% dist = Distancia entre los puntos del mallado regular interior para el
% cálculo de la solución particular u_P en cada problema elíptico. Se usará
% al invocar la función MSF_MSPEliptDirGen.
% Nb = Número de puntos de la frontera para aplicar el MSF en cada problema
% elíptico. Se usará al invocar la función MSF_MSPEliptDirGen.
% lambda0 = Factor para ubicar las fuentes fuera de la bola. Se usará al 
% invocar la función MSF_MSPEliptDirGen.
% lambda1 = Constante que aparece en la definición de las funciones
% radiales f y F, relacionadas por Lap(F) = f. Se usará al invocar la función 
% MSF_MSPEliptDirGen.
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

% Discretización en tiempo.

t = linspace(0,T,M+1);
k = T/M;
k1 = 1/k;
a = @(x) k1;

% Inicialización.

u = cell(1,M+1);
u{1} = u0;

% Evolución temporal.

for n = 1:M
    un = u{n};
    tn1 = t(n+1);
    h1 = @(x) h(x,tn1);
    h_elipt = @(x) h1(x) + k1*un(x);
    g1 = @(x) g(x,tn1);
    u{n+1} = MSF_MSPEliptDirGen(x1_eq,x2_eq,xc,dist,Nb,lambda0,lambda1,a,h_elipt,g1);
end

% Visualización.

% Generamos una malla de puntos dentro del dominio.

theta_plot = linspace(0,2*pi,300);
x1_b = x1_eq(theta_plot);
x2_b = x2_eq(theta_plot);
x1_min = min(x1_b); x1_max = max(x1_b);
x2_min = min(x2_b); x2_max = max(x2_b);

Nmalla = 300;
[X,Y] = meshgrid(linspace(x1_min,x1_max,Nmalla),linspace(x2_min,x2_max,Nmalla));
in = inpolygon(X,Y,x1_b,x2_b);

% Graficamos las soluciones y el error para cada instante.

for n = 1:M+1
    U = nan(size(X));
    Ue = nan(size(X));
    E = nan(size(X));
    for i = 1:size(X,1)
        for j = 1:size(X,2)
            if in(i,j)
                x = [X(i,j); Y(i,j)];
                U(i,j) = u{n}(x);
                Ue(i,j) = u_exacta(x,t(n));
                E(i,j) = abs(U(i,j) - Ue(i,j));
            end
        end
    end

    % Solución numérica

    figure(1)
    surf(X,Y,U,'EdgeColor','none');
    colorbar;
    title(['Solución del método en t = ', num2str(t(n))]);
    xlabel('x'); ylabel('y'); zlabel('u(x)');
    axis equal tight;
    view(3);
    shading interp;
    pause(0.1);

    % Solución exacta
    figure(2)
    surf(X,Y,Ue,'EdgeColor','none');
    colorbar;
    title(['Solución exacta en t = ', num2str(t(n))]);
    xlabel('x'); ylabel('y'); zlabel('u_{exacta}(x)');
    axis equal tight;
    view(3);
    shading interp;
    pause(0.1);

    % Error absoluto
    figure(3)
    surf(X,Y,E,'EdgeColor','none');
    colorbar;
    title(['Error absoluto en t = ', num2str(t(n))]);
    xlabel('x'); ylabel('y'); zlabel('|u - u_{exacta}|');
    axis equal tight;
    view(3);
    shading interp;
    pause(0.1);
end

end
