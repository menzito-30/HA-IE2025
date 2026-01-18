%Optimized for speed by SeHyoun Ahn
% ---
% TÍTULO: huggett_equilibrium_iterate.m (basado en la estructura)
% AUTOR: Achdou, Han, Lasry, Lions, Moll (con optimizaciones de SeHyoun Ahn)
% PROPÓSITO: Resuelve el modelo de Aiyagari-Bewley-Huggett en tiempo continuo
%            y ENCUENTRA la tasa de interés de equilibrio 'r'
%            que vacía el mercado de activos (S(r) = 0).
% MÉTODO: Bisección (Bisection Search) para el equilibrio 'r'.
%         Dentro de cada iteración de 'r', usa:
%         1. Iteración de Función de Valor (VFI) con Método Implícito para la HJB.
%         2. Solución de sistema lineal (A'g=0) para la Kolmogorov Forward (KF).
% ---

clear all; clc; close all;

tic; % Inicia el cronómetro

% ---
% 1. PARÁMETROS ECONÓMICOS
% ---
s = 2; % Coeficiente de aversión al riesgo relativo (γ)
rho = 0.05; % Tasa de descuento del tiempo (ρ)
z1 = .1; % Estado de bajos ingresos
z2 = .2; % Estado de altos ingresos
z = [z1,z2]; % Vector de estados de ingreso
la1 = 1.2; % Tasa de transición z1 -> z2 (λ_1)
la2 = 1.2; % Tasa de transición z2 -> z1 (λ_2)
la = [la1,la2]; % Vector de tasas de transición

% ---
% 2. PARÁMETROS DE BÚSQUEDA DE EQUILIBRIO (BISECCIÓN)
% ---
r0 = 0.03; % Guess inicial para la tasa de interés 'r'
rmin = 0.01; % Límite inferior para la búsqueda de 'r'
rmax = 0.04; % Límite superior para la búsqueda de 'r'
Ir = 40; % Máximo de iteraciones para la búsqueda de 'r' (bisección)
crit_S = 10^(-5); % Criterio de convergencia para el mercado (S(r) ≈ 0)

% ---
% 3. PARÁMETROS NUMÉRICOS (GRILLA Y VFI)
% ---
I= 1000; % Número de puntos en la grilla de activos
amin = -0.15; % Límite de endeudamiento (a_barra)
amax = 5; % Límite superior de la grilla de activos
a = linspace(amin,amax,I)'; % Vector (columna) de la grilla de activos
da = (amax-amin)/(I-1); % Tamaño del paso en la grilla

% Matrices auxiliares para operaciones vectorizadas
aa = [a,a]; % Repite la grilla de activos para ambos estados z
zz = ones(I,1)*z; % Repite los estados de ingreso para toda la grilla 'a'

% Parámetros de la VFI (Método Implícito)
maxit= 100; % Máximo de iteraciones para la HJB (VFI)
crit = 10^(-6); % Criterio de convergencia para la HJB
Delta = 1000; % Paso de tiempo (grande) para el método implícito
              % (un Delta grande acelera la convergencia)

% Pre-asignación de matrices
dVf = zeros(I,2); % Derivada hacia adelante (forward)
dVb = zeros(I,2); % Derivada hacia atrás (backward)
c = zeros(I,2); % Función de consumo

% Matriz (sparse) que captura los saltos entre estados z_1 y z_2
Aswitch = [-speye(I)*la(1),speye(I)*la(1);speye(I)*la(2),-speye(I)*la(2)];

% ---
% 4. GUESS INICIAL PARA LA FUNCIÓN DE VALOR
% ---
r = r0; % Empezar la búsqueda con el guess inicial
% Guess v0: valor de consumir la dotación + intereses para siempre
v0(:,1) = (z(1) + r.*a).^(1-s)/(1-s)/rho;
v0(:,2) = (z(2) + r.*a).^(1-s)/(1-s)/rho;

% ---
% 5. LOOP EXTERNO: BÚSQUEDA DE EQUILIBRIO 'r' (BISECCIÓN)
% ---
% Este loop itera para encontrar la 'r' que vacía el mercado (S(r)=0)
for ir=1:Ir;

% Guardar historial de la búsqueda de 'r' (opcional)
r_r(ir)=r;
rmin_r(ir)=rmin;
rmax_r(ir)=rmax;
    
% "Warm Start": Usar la solución de V de la iteración anterior
% (para la 'r' anterior) como guess inicial. Acelera la convergencia.
if ir>1
v0 = V_r(:,:,ir-1);
end

v = v0; % Inicializar la VFI con el guess

% ---
% 6. LOOP INTERNO: RESOLVER LA HJB (VFI) para una 'r' DADA
% ---
% Este loop resuelve el problema del hogar para una 'r' fija
for n=1:maxit
    V = v;
    V_n(:,:,n)=V; % Guardar historial de VFI (opcional)
    
    % ---
    % 6A. CÁLCULO DE DERIVADAS (Esquema "Upwind")
    % ---
    
    % Derivada hacia adelante
    dVf(1:I-1,:) = (V(2:I,:)-V(1:I-1,:))/da;
    dVf(I,:) = (z + r.*amax).^(-s); % Condición de contorno en amax
    
    % Derivada hacia atrás
    dVb(2:I,:) = (V(2:I,:)-V(1:I-1,:))/da;
    dVb(1,:) = (z + r.*amin).^(-s); % Condición de contorno en amin
    
    I_concave = dVb > dVf; % Indicador de concavidad
    
    % Calcular consumo y ahorro (drift) para ambas derivadas
    % Se usa max(..., 10^-10) para estabilidad numérica (evitar dVf<=0)
    cf = max(dVf,10^(-10)).^(-1/s); % Consumo con v_f'
    ssf = zz + r.*aa - cf; % Ahorro (drift) con v_f'
    cb = max(dVb,10^(-10)).^(-1/s); % Consumo con v_b'
    ssb = zz + r.*aa - cb; % Ahorro (drift) con v_b'
    
    % Consumo si el drift es cero (steady state individual)
    c0 = zz + r.*aa;
    
    % SELECCIÓN UPWIND
    % Elegir la derivada basada en la dirección del ahorro (drift)
    If = ssf > 0; % Usar forward diff si el ahorro (drift) es > 0
    Ib = ssb < 0; % Usar backward diff si el ahorro (drift) es < 0
    I0 = (1-If-Ib); % Usar c0 si el ahorro (drift) es = 0
    
    % Esta es una implementación alternativa del upwind:
    % selecciona el *consumo* (cf, cb, c0) en lugar de la *derivada* (dVf, dVb, dV0)
    c = cf.*If + cb.*Ib + c0.*I0;
    u = c.^(1-s)/(1-s); % Utilidad
    
    % ---
    % 6B. CONSTRUCCIÓN DE LA MATRIZ GENERADORA 'A'
    % ---
    
    % Componentes de la matriz tridiagonal (del drift s(a)*v'(a))
    X = -min(ssb,0)/da;
    Y = -max(ssf,0)/da + min(ssb,0)/da;
    Z = max(ssf,0)/da;
    
    % Ensamblar A para cada estado 'z'
    A1=spdiags(Y(:,1),0,I,I)+spdiags(X(2:I,1),-1,I,I)+spdiags([0;Z(1:I-1,1)],1,I,I);
    A2=spdiags(Y(:,2),0,I,I)+spdiags(X(2:I,2),-1,I,I)+spdiags([0;Z(1:I-1,2)],1,I,I);
    % Combinar con los saltos (Aswitch) para la matriz A (2I x 2I)
    % A = Operador HJB (drift + saltos de Poisson)
    A = [A1,sparse(I,I);sparse(I,I),A2] + Aswitch;
    
    % Chequeo: Las filas de la matriz A deben sumar 0
    if max(abs(sum(A,2)))>10^(-9)
       disp('Improper Transition Matrix')
       break
    end
    
    % ---
    % 6C. RESOLVER LA HJB (Método Implícito)
    % ---
    
    % La iteración es: ( (1/Δ + ρ)I - A ) V_new = u + V_old/Δ
    % B = ( (1/Δ + ρ)I - A )
    B = (1/Delta + rho)*speye(2*I) - A;

    % Apilar vectores u y V (de la iteración n)
    u_stacked = [u(:,1);u(:,2)];
    V_stacked = [V(:,1);V(:,2)];
    
    % b = u + V_old/Δ
    b = u_stacked + V_stacked/Delta;
    % Resolver B * V_new = b  (este es el paso VFI)
    V_stacked = B\b; %SOLVE SYSTEM OF EQUATIONS
    
    V = [V_stacked(1:I),V_stacked(I+1:2*I)]; % Des-apilar V_new
    
    % Chequeo de convergencia de la VFI
    Vchange = V - v;
    v = V; % Actualizar v para la siguiente iteración

    dist(n) = max(max(abs(Vchange)));
    if dist(n)<crit
        disp('Value Function Converged, Iteration = ')
        disp(n)
        break
    end
end % Fin del loop VFI (interno)
toc; % Reportar tiempo de VFI para esta 'r'

%%%%%%%%%%%%%%%%%%%%%%%%%%
% FOKKER-PLANCK EQUATION %
%%%%%%%%%%%%%%%%%%%%%%%%%%

% ---
% 7. RESOLVER LA KOLMOGOROV FORWARD (KF)
% ---
% Una vez que V converge, 'A' es el generador final.
% La distribución estacionaria 'g' resuelve A'g = 0.
AT = A'; % Usar la transpuesta de la matriz A (el operador adjunto)
b = zeros(2*I,1); % Buscamos el vector en el kernel (null space)

% "Dirty fix" para la singularidad
% A'g=0 tiene infinitas soluciones (g=0, g=c*g_ss).
% Fijamos un punto (g(1)=0.1) y luego normalizamos la suma total a 1.
i_fix = 1;
b(i_fix)=.1;
row = [zeros(1,i_fix-1),1,zeros(1,2*I-i_fix)];
AT(i_fix,:) = row;

%Solve linear system
gg = AT\b; % Solución 'g' no normalizada
g_sum = gg'*ones(2*I,1)*da; % Calcular la integral total (suma de g_i * da)
gg = gg./g_sum; % Normalizar para que la integral sea 1

g = [gg(1:I),gg(I+1:2*I)]; % Des-apilar la distribución 'g'

% Chequeos (opcionales)
check1 = g(:,1)'*ones(I,1)*da; % Masa total tipo 1
check2 = g(:,2)'*ones(I,1)*da; % Masa total tipo 2

% Guardar resultados para esta iteración de 'r'
g_r(:,:,ir) = g; % Guardar distribución
adot(:,:,ir) = zz + r.*aa - c; % Guardar política de ahorro
V_r(:,:,ir) = V; % Guardar función de valor

% ---
% 8. CHEQUEO DE EQUILIBRIO Y ACTUALIZACIÓN (BISECCIÓN)
% ---
% Calcular la oferta agregada de activos S(r) = E[a]
S(ir) = g(:,1)'*a*da + g(:,2)'*a*da;

% Algoritmo de Bisección
if S(ir)>crit_S % Exceso de Oferta (S > 0)
    disp('Excess Supply')
    rmax = r; % El 'r' actual es demasiado alto, se vuelve el nuevo máximo
    r = 0.5*(r+rmin); % Probar un 'r' a medio camino entre el actual y el mínimo
elseif S(ir)<-crit_S; % Exceso de Demanda (S < 0)
    disp('Excess Demand')
    rmin = r; % El 'r' actual es demasiado bajo, se vuelve el nuevo mínimo
    r = 0.5*(r+rmax); % Probar un 'r' a medio camino entre el actual y el máximo
elseif abs(S(ir))<crit_S; % Equilibrio Encontrado (S ≈ 0)
    display('Equilibrium Found, Interest rate =')
    disp(r)
    break % Salir del loop de bisección
end

end % Fin del loop de bisección (externo)

% ---
% 9. GRÁFICOS (para la 'r' de equilibrio encontrada)
% ---
amax1 = 0.6; % Límite superior del eje x para el zoom
amin1 = amin-0.03; % Límite inferior del eje x para el zoom

% Gráfico 1: Política de Ahorro
figure(1)
set(gca,'FontSize',16)
h1 = plot(a,adot(:,1,ir),'b',a,adot(:,2,ir),'r',linspace(amin1,amax1,I),zeros(1,I),'k--','LineWidth',2)
legend(h1,'s_1(a)','s_2(a)','Location','NorthEast')
text(-0.155,-.105,'$\underline{a}$','FontSize',16,'interpreter','latex') % Etiqueta para a_min
line([amin amin], [-.1 .08],'Color','Black','LineStyle','--') % Línea vertical en a_min
xlabel('Wealth, $a$','interpreter','latex')
ylabel('Savings, $s_i(a)$','interpreter','latex')
xlim([amin1 amax1])

% Gráfico 2: Distribución de Riqueza
figure(2)
set(gca,'FontSize',16)
h1 = plot(a,g_r(:,1,ir),'b',a,g_r(:,2,ir),'r','LineWidth',2)
legend(h1,'g_1(a)','g_2(a)')
text(-0.155,-.12,'$\underline{a}$','FontSize',16,'interpreter','latex') % Etiqueta para a_min
line([amin amin], [0 max(max(g_r(:,:,ir)))],'Color','Black','LineStyle','--') % Línea vertical en a_min
xlabel('Wealth, $a$','interpreter','latex')
ylabel('Densities, $g_i(a)$','interpreter','latex')
xlim([amin1 amax1])
ylim([0 3])