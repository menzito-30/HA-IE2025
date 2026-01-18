clear all; clc; close all;

tic; % Inicia el cronómetro

% ---
% 1. PARÁMETROS DEL MODELO
% ---

% Parámetro de preferencias (CRRA)
% NOTA: 's' aquí es gamma (γ) en la notación estándar del paper,
% el coeficiente de aversión al riesgo relativo.
s = 2;
rho = 0.05; % Tasa de descuento del tiempo (ρ)

% Proceso de ingreso (estocástico)
z1 = .1; % Estado de bajos ingresos (ej. desempleado)
z2 = .2; % Estado de altos ingresos (ej. empleado)
z = [z1,z2]; % Vector de estados de ingreso

% Tasas de transición de Poisson (λ)
% la1 es la tasa de z1 -> z2 (ej. encontrar trabajo)
% la2 es la tasa de z2 -> z1 (ej. perder trabajo)
% (Los valores originales están comentados, se usan los de abajo)
% la1 = 0.001;
% la2 = 0.001;
% la1 = 0.02;
% la2 = 0.03;
la1 = 1.5;
la2 = 1;
la = [la1,la2];


% ---
% 2. PARÁMETROS NUMÉRICOS Y GRILLA
% ---

% Grilla de activos (a)
I=500; % Número de puntos en la grilla de activos
amin = -0.15; % Límite de endeudamiento (restricción, a_barra)
amax = 10; % Límite superior (truncamiento) de la grilla de activos
a = linspace(amin,amax,I)'; % Vector columna de la grilla de activos
da = (amax-amin)/(I-1); % Tamaño del paso (step size)

% Matrices auxiliares para 'broadcast'
aa = [a,a]; % Matriz de grilla de activos (I x 2)
zz = ones(I,1)*z; % Matriz de estados de ingreso (I x 2)


% Parámetros del iterador de la HJB (Método Implícito)
maxit= 100; % Máximo de iteraciones para la HJB
crit = 10^(-6); % Criterio de convergencia (tolerancia)
Delta = 40; % Paso de tiempo para el método implícito [ver Apéndice Numérico, ec. 16,  193]
           % Un Delta grande acelera la convergencia (similar a T->infinito)

% Pre-asignación de matrices (buena práctica en MATLAB)
dVf = zeros(I,2); % Derivada hacia adelante (forward) de V
dVb = zeros(I,2); % Derivada hacia atrás (backward) de V
c = zeros(I,2); % Consumo

% Matriz de transición de estado de ingreso (saltos de Poisson)
% Esta es la parte "off-diagonal" del operador A
% Representa los términos λ_j(v_{-j} - v_j)
Aswitch = [-speye(I)*la(1),speye(I)*la(1);speye(I)*la(2),-speye(I)*la(2)];

% Grilla de la tasa de interés (r)
% El script resolverá el modelo para cada una de estas 'r'
Ir = 20; % Número de puntos en la grilla de 'r'
% rmin = 0.00001;
rmin = -0.05;
rmax = 0.04;
r_grid = linspace(rmin,rmax,Ir);

% ---
% 3. VALOR INICIAL (GUESS) PARA LA FUNCIÓN DE VALOR
% ---
r = r_grid(1);
% El valor inicial es el valor de consumir la dotación + intereses para siempre
% (Se usa max(r, 0.01) para evitar problemas si r es negativo)
v0(:,1) = (z(1) + max(r,0.01).*a).^(1-s)/(1-s)/rho;
v0(:,2) = (z(2) + max(r,0.01).*a).^(1-s)/(1-s)/rho;

% ---
% 4. LOOP PRINCIPAL: RESOLVER PARA CADA TASA DE INTERÉS 'r'
% ---

for ir=1:Ir;

r = r_grid(ir);

% "Warm start": Usar la solución del 'r' anterior como guess inicial
% Acelera mucho la convergencia.
if ir>1
v0 = V_r(:,:,ir-1);
end

v = v0; % Inicializar la función de valor para este loop de 'r'

% ---
% 4A. LOOP INTERNO: RESOLVER LA HJB (Iteración de Función de Valor)
% ---
for n=1:maxit
    V = v; % Guardar el valor de la iteración anterior
    V_n(:,:,n)=V; % (Opcional: guardar el historial de iteraciones)
    
    % ---
    % CÁLCULO DE DERIVADAS (Esquema "Upwind")
    % ---
    
    % 1. Derivada hacia adelante (Forward difference)
    dVf(1:I-1,:) = (V(2:I,:)-V(1:I-1,:))/da;
    dVf(I,:) = (z + r.*amax).^(-s); % Condición de contorno en amax (state constraint)
    
    % 2. Derivada hacia atrás (Backward difference)
    dVb(2:I,:) = (V(2:I,:)-V(1:I-1,:))/da;
    % Condición de contorno en amin (restricción de endeudamiento)
    % v'(amin) = u'(z + r*amin) [Apéndice Numérico, ec. 7, cite: 23]
    dVb(1,:) = (z + r.*amin).^(-s); %state constraint boundary condition
    
    I_concave = dVb > dVf; %indicator whether value function is concave (problems arise if this is not the case)
    
    % 3. Calcular consumo y ahorro (drift) para ambas derivadas
    % Forward
    cf = dVf.^(-1/s);
    ssf = zz + r.*aa - cf; % Ahorro (drift s(a)) si se usa V_f
    % Backward
    cb = dVb.^(-1/s);
    ssb = zz + r.*aa - cb; % Ahorro (drift s(a)) si se usa V_b
    %consumption and derivative of value function at steady state
    % 4. Consumo en el "steady state" (drift = 0)
    c0 = zz + r.*aa;
    dV0 = c0.^(-s);
    
    % dV_upwind makes a choice of forward or backward differences based on
    % the sign of the drift    
    
    % 5. SELECCIÓN UPWIND
    % El esquema "upwind" elige la derivada basada en la dirección
    % del ahorro (el "drift" de la variable de estado 'a').
    If = ssf > 0; % Usar forward diff si el drift (ahorro) es positivo
    Ib = ssb < 0; % Usar backward diff si el drift (ahorro) es negativo
    I0 = (1-If-Ib); % Usar V_0 si el drift es cero
    %make sure backward difference is used at amax
    %Ib(I,:) = 1; If(I,:) = 0;
    %STATE CONSTRAINT at amin: USE BOUNDARY CONDITION UNLESS muf > 0:
    %already taken care of automatically
    
    % 6. Combinar para obtener la derivada "upwind" correcta
    dV_Upwind = dVf.*If + dVb.*Ib + dV0.*I0; %important to include third term
    c = dV_Upwind.^(-1/s); % Política de consumo final
    u = c.^(1-s)/(1-s); % Utilidad
    
    % ---
    % 4B. CONSTRUCCIÓN DE LA MATRIZ GENERADORA 'A' (Operador HJB)
    % ---
    
    % Esta matriz 'A' es el generador infinitesimal del proceso estocástico.
    % Es una matriz sparse (dispersa) que representa el operador HJB
    % discretizado. [Apéndice Numérico, ec. 15, cite: 127]
    
    % Componentes de la matriz tridiagonal (para la parte del drift s(a)*V')
    X = - min(ssb,0)/da;
    Y = - max(ssf,0)/da + min(ssb,0)/da;
    Z = max(ssf,0)/da;
    
    % Ensamblar las matrices tridiagonales para cada estado de ingreso
    A1=spdiags(Y(:,1),0,I,I)+spdiags(X(2:I,1),-1,I,I)+spdiags([0;Z(1:I-1,1)],1,I,I);
    A2=spdiags(Y(:,2),0,I,I)+spdiags(X(2:I,2),-1,I,I)+spdiags([0;Z(1:I-1,2)],1,I,I);
    
    % Ensamblar la matriz A completa (2I x 2I)
    % Combina el drift (A1, A2) con los saltos de ingreso (Aswitch)
    A = [A1,sparse(I,I);sparse(I,I),A2] + Aswitch;
    
    % Chequeo de sanidad: Las filas de un generador 'A' deben sumar 0
    if max(abs(sum(A,2)))>10^(-9)
       disp('Improper Transition Matrix')
       break
    end
    
    % ---
    % 4C. RESOLVER LA ECUACIÓN DE BELLMAN (Método Implícito)
    % ---
    
    % La iteración es: ( (1/Δ + ρ)I - A ) V^{n+1} = u^n + V^n/Δ
    % Sea B = ( (1/Δ + ρ)I - A )
    % Sea b = u^n + V^n/Δ
    % Resolvemos el sistema lineal B * V^{n+1} = b
    % [Ver Apéndice Numérico, ec. 16, cite: 193]
    
    B = (1/Delta + rho)*speye(2*I) - A;
    
    % Apilar vectores u y V para el sistema (2I x 1)
    u_stacked = [u(:,1);u(:,2)];
    V_stacked = [V(:,1);V(:,2)];
    
    b = u_stacked + V_stacked/Delta;
    V_stacked = B\b; %SOLVE SYSTEM OF EQUATIONS (paso clave de VFI)
    
    V = [V_stacked(1:I),V_stacked(I+1:2*I)]; % Des-apilar V^{n+1}
    
    Vchange = V - v;
    v = V; % Actualizar la función de valor para la próxima iteración

    dist(n) = max(max(abs(Vchange)));
    if dist(n)<crit
        disp('Value Function Converged, Iteration = ')
        disp(n)
        break
    end
end
toc; % Muestra el tiempo para la convergencia de la HJB

%%%%%%%%%%%%%%%%%%%%%%%%%%
% FOKKER-PLANCK EQUATION %
%%%%%%%%%%%%%%%%%%%%%%%%%%

% ---
% 5. RESOLVER LA ECUACIÓN DE KOLMOGOROV FORWARD (KF)
% ---

% Una vez que HJB converge, 'A' es el generador de estado estacionario final.
% La distribución estacionaria 'g' resuelve A'g = 0.
% Esta es la propiedad clave de "transposición" (adjoint). [cite: 245, 247]

AT = A'; % Transponer la matriz generadora
b = zeros(2*I,1); % Vector de ceros (buscamos A'g = 0)

%need to fix one value, otherwise matrix is singular
% "Dirty fix" para la singularidad [Apéndice Numérico, p. 9, cite: 251, 252, 253]
% A' es singular (las filas suman 0), por lo que A'g=0 tiene infinitas
% soluciones. Fijamos un punto de la distribución (ej. g(1)=0.1)
% para anclar la solución, y luego normalizamos.
i_fix = 1;
b(i_fix)=.1;
row = [zeros(1,i_fix-1),1,zeros(1,2*I-i_fix)];
AT(i_fix,:) = row;

%Solve linear system
gg = AT\b; % Resolver el sistema lineal para 'g' (no normalizado)
g_sum = gg'*ones(2*I,1)*da; % Integral = suma(g_i * da)
gg = gg./g_sum; % Normalizar la distribución para que sume (integre) a 1

g = [gg(1:I),gg(I+1:2*I)]; % Des-apilar el vector de distribución 'g'

% (Opcional) Chequear que las densidades marginales suman correctamente
check1 = g(:,1)'*ones(I,1)*da;
check2 = g(:,2)'*ones(I,1)*da;

% ---
% 6. GUARDAR RESULTADOS Y CONTINUAR EL LOOP DE 'r'
% ---
g_r(:,:,ir) = g; % Guardar la distribución para este 'r'
adot(:,:,ir) = zz + r.*aa - c; % Guardar la política de ahorro (drift) para este 'r'
V_r(:,:,ir) = V; % Guardar la función de valor para este 'r'
dV_r(:,:,ir) = dV_Upwind; % Guardar la derivada de V (política de consumo) para este 'r'
c_r(:,:,ir) = c; % Guardar la política de consumo para este 'r'

% CALCULAR LA OFERTA AGREGADA DE ACTIVOS S(r)
% S(r) = E[a] = integral(a * g(a) da) [cite: 280]
% Esta es la métrica clave que produce este script.
S(ir) = g(:,1)'*a*da + g(:,2)'*a*da;
end % Fin del loop principal (sobre 'r')

% ---
% ---
% ---
% ---
% 7. GRÁFICOS 
% ---
% Esta sección genera los gráficos que se encuentran en el Apéndice Numérico.
%
% NOTA: Los índices 'ir' están fijados (hardcoded) en el script original
% para replicar los gráficos del paper.
% ---

% ir = 9;
% set(gca,'FontSize',14)
% h1 = plot(a,dV_r(:,:,ir),a,dV_r(:,:,ir+1),'LineWidth',2)
% legend(h1,'v_1\prime(a,r1)','v_2(a,r1)','v_1(a,r2)','v_2(a,r2)')
% grid
% xlabel('a')
% ylabel('s_i(a)')
% xlim([amin amax])

% ---
% Gráfico 1: Función de Política de Ahorro s_j(a)
% Compara la política de ahorro para dos tasas de interés (ir y ir+1)
% ---
%{     
    ir = 1; % Fija el índice de la primera tasa de interés para el gráfico
    set(gca,'FontSize',14) % Fija el tamaño de la fuente
    % Dibuja 5 líneas: 2 para adot(ir), 2 para adot(ir+1) y 1 para la línea cero
    h1 = plot(a,adot(:,:,ir),a,adot(:,:,ir+1),a,zeros(1,I),'--','LineWidth',2);
    % NOTA: El siguiente comando 'legend' genera el "Warning" que viste.
    % Solo proporciona 4 etiquetas para las 5 líneas dibujadas, ignorando la línea 'zeros'.
    legend(h1,'s_1(a,r1)','s_2(a,r1)','s_1(a,r2)','s_2(a,r2)')
    grid
    xlabel('a')
    ylabel('s_i(a)')
    xlim([amin amax])
%}

% ---
% Gráfico 2: Función de Política de Consumo c_j(a)
% ---
%{    
    ir = 9; % Fija el índice de la tasa de interés para este gráfico
    set(gca,'FontSize',14)
    % Dibuja 5 líneas: 2 para c_r(ir), 2 para c_r(ir+1) y 1 para la línea cero
    h1 = plot(a,c_r(:,:,ir),a,c_r(:,:,ir+1),a,zeros(1,I),'--','LineWidth',2);
    % NOTA: Igual que antes, este 'legend' ignora la línea 'zeros', causando el Warning.
    legend(h1,'c_1(a,r1)','c_2(a,r1)','c_1(a,r2)','c_2(a,r2)')
    grid
    xlabel('a')
    ylabel('c_i(a)')
    xlim([amin amax])
%}
    
% ---
% Gráfico 3: Densidades de Riqueza Estacionarias g_j(a)
% ---
%{
    figure; % Crea una nueva figura
    amax1 = .8; % Fija un límite superior del eje x (para hacer zoom)
    set(gca,'FontSize',14)
    % Dibuja las 2 densidades (g1 y g2) para la tasa 'ir' fijada en la línea 193
    h1 = plot(a,g_r(:,:,ir),'LineWidth',2);
    legend(h1,'g_1(a)','g_2(a)')
    grid
    xlabel('a')
    ylabel('g_i(a)')
    xlim([amin amax1]) % Aplica el zoom
%}


% ---
% Gráfico 4: Curva de Oferta de Activos S(r)
% ---
figure; % Crea una nueva figura

% Prepara los ejes para las líneas de referencia
Smax = max(S);
amin1 = 1.1*amin;
aaa = linspace(amin1,Smax,Ir); % Eje x para la línea r=rho
rrr = linspace(rmin,0.06,Ir); % Eje y para las líneas S=0 y a=a_barra

set(gca,'FontSize',14)
% Dibuja la curva S(r) y las 3 líneas de referencia
plot(S,r_grid,zeros(1,Ir),rrr,zeros(1,Ir)+amin,rrr,'--',aaa,ones(1,Ir)*rho,'--','LineWidth',2)

% Añade etiquetas de texto manualmente usando coordenadas
text(-0.1,0.045,'$r = \rho$','FontSize',16,'interpreter','latex')
text(-0.05,0,'$S(r)$','FontSize',16,'interpreter','latex')
text(-0.145,0.01,'$a=\underline{a}$','FontSize',16,'interpreter','latex')

% Etiquetas de los ejes
ylabel('$r$','FontSize',16,'interpreter','latex')
xlabel('$S(r)$','FontSize',16,'interpreter','latex')

% Límites de los ejes
ylim([rmin 0.06])
xlim([amin1 Smax])
% print -depsc asset_supply.eps % (Comentado para no guardar archivo)