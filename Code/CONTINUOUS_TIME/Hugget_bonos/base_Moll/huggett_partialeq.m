%Optimized for speed by SeHyoun Ahn

clear all; clc; close all;

% Inicia el cronómetro para medir el tiempo de ejecución
tic;

% -------------------------------------------------------------------------
% 1. DEFINICIÓN DE PARÁMETROS DEL MODELO
% -------------------------------------------------------------------------

s = 1.2;       % Coeficiente de aversión relativa al riesgo (gamma en los apuntes) 
r = 0.035;     % Tasa de interés (fija, estamos en equilibrio parcial) 
rho = 0.05;    % Tasa de descuento subjetiva (impaciencia) 
z1 = .1;       % Estado de ingresos bajos 
z2 = .2;       % Estado de ingresos altos 
z = [z1,z2];   % Vector de estados de ingreso
la1 = 1.5;     % Intensidad de Poisson para pasar de z1 a z2 (lambda_1) 
la2 = 1;       % Intensidad de Poisson para pasar de z2 a z1 (lambda_2) 
la = [la1,la2]; % Vector de intensidades

% -------------------------------------------------------------------------
% 2. CREACIÓN DEL GRID (REJILLA)
% -------------------------------------------------------------------------
% El método de diferencias finitas requiere discretizar el espacio de estados 

I=500;             % Número de puntos en la grilla de activos 'a'
amin = -0.02;      % Límite inferior de activos (restricción de endeudamiento, \underline{a}) 
amax = 3;          % Límite superior de activos (truncamiento artificial)
a = linspace(amin,amax,I)'; % Vector columna de la grilla de activos (equiespaciada)
da = (amax-amin)/(I-1);  % Distancia entre puntos de la grilla (\Delta a) 

% Matrices auxiliares para cálculos vectorizados
aa = [a,a];        % Repite la grilla de activos para ambos estados de 'z'
zz = ones(I,1)*z;  % Matriz (I x 2) con los valores de 'z' correspondientes

% -------------------------------------------------------------------------
% 3. PARÁMETROS DE ITERACIÓN NUMÉRICA
% -------------------------------------------------------------------------

maxit= 100;        % Número máximo de iteraciones para la función de valor
crit = 10^(-6);    % Criterio de convergencia (tolerancia)
Delta = 1000;      % Tamaño del paso para el método implícito 
                   % Un Delta grande se acerca a la solución estacionaria rho*v = ... 

% Inicialización de matrices para las derivadas y la política de consumo
dVf = zeros(I,2); % Derivada aproximada hacia adelante (Forward difference)
dVb = zeros(I,2); % Derivada aproximada hacia atrás (Backward difference)
c = zeros(I,2);   % Política de consumo

% Matriz (sparse) que captura los saltos entre estados de ingreso (z1 <-> z2)
% Esto representa los términos \lambda_j(v_{-j} - v_j) en la HJB 
% Es una matriz de bloque (2I x 2I) 
Aswitch = [-speye(I)*la(1),speye(I)*la(1);speye(I)*la(2),-speye(I)*la(2)];

% -------------------------------------------------------------------------
% 4. CONJETURA INICIAL (INITIAL GUESS) PARA LA FUNCIÓN DE VALOR
% -------------------------------------------------------------------------
% Una conjetura común es el valor de consumir permanentemente y+ra 
v0(:,1) = (z(1) + r.*a).^(1-s)/(1-s)/rho;
v0(:,2) = (z(2) + r.*a).^(1-s)/(1-s)/rho;

v = v0; % 'v' es la función de valor que se actualizará en cada iteración

% -------------------------------------------------------------------------
% 5. BUCLE PRINCIPAL: ITERACIÓN DE LA FUNCIÓN DE VALOR (HJB)
% -------------------------------------------------------------------------
% Se itera hasta que la función de valor converja 

for n=1:maxit
    V = v; % Guarda la función de valor de la iteración anterior (v^n)
    V_n(:,:,n)=V; % (Opcional: guarda el historial de iteraciones)
    
    % --- 5.1. Aproximación de derivadas (Diferencias Finitas) ---
    
    % Derivada hacia adelante: v'(a_i) \approx (v_{i+1} - v_i) / da 
    dVf(1:I-1,:) = (V(2:I,:)-V(1:I-1,:))/da;
    dVf(I,:) = (z + r.*amax).^(-s); % Condición en el límite superior (amax)
    
    % Derivada hacia atrás: v'(a_i) \approx (v_i - v_{i-1}) / da 
    dVb(2:I,:) = (V(2:I,:)-V(1:I-1,:))/da;
    dVb(1,:) = (z + r.*amin).^(-s); % Condición de restricción estatal en amin 
                                   % v'(\underline{a}) = u'(z + r*\underline{a})
    
    % Comprobación de concavidad (útil para depurar) 
    I_concave = dVb > dVf; 
    
    % --- 5.2. ESQUEMA "UPWIND" ---
    % El esquema "upwind" elige la diferencia (adelante/atrás) basándose
    % en la dirección del "drift" (ahorro, s_j(a)) 
    
    % Primero, calcula el drift (ahorro) para AMBAS aproximaciones 
    
    % 1. Ahorro con derivada hacia adelante (forward)
    cf = dVf.^(-1/s);      % c = (u')^{-1}(v') => c = (v')^{-1/s} 
    ssf = zz + r.*aa - cf; % s_j(a) = z_j + ra - c_j(a) (drift) 
    
    % 2. Ahorro con derivada hacia atrás (backward)
    cb = dVb.^(-1/s);
    ssb = zz + r.*aa - cb;
    
    % 3. Ahorro nulo (c_j(a) = z_j + ra)
    % Se usa si el drift cambia de signo (s_F <= 0 <= s_B) 
    c0 = zz + r.*aa;
    dV0 = c0.^(-s);
    
    % Indicadores para elegir la derivada correcta 
    If = ssf > 0; % Usar 'forward' si el drift (ahorro) es positivo
    Ib = ssb < 0; % Usar 'backward' si el drift (ahorro) es negativo
    I0 = (1-If-Ib); % Usar 'cero' (dV0) si el drift está entre ssf y ssb
    
    % Selecciona la derivada apropiada según el esquema upwind
    dV_Upwind = dVf.*If + dVb.*Ib + dV0.*I0; 
    
    % Calcula la política de consumo final y la utilidad 
    c = dV_Upwind.^(-1/s);
    u = c.^(1-s)/(1-s);
    
    % --- 5.3. CONSTRUCCIÓN DE LA MATRIZ DE TRANSICIÓN 'A' ---
    % Esta matriz 'A' es el "generador infinitesimal" discretizado 
    % Representa la parte de la derivada espacial (v') de la HJB.
    % Los términos vienen de reordenar la HJB discretizada 
    
    % x_ij: coeficiente para v_{i-1,j}
    X = - min(ssb,0)/da;
    % y_ij: coeficiente para v_{i,j}
    Y = - max(ssf,0)/da + min(ssb,0)/da;
    % z_ij: coeficiente para v_{i+1,j}
    Z = max(ssf,0)/da;
    
    % (Nota: los términos -lambda_j no están en Y
    % porque se añaden por separado a través de 'Aswitch')

    % Construye bloques diagonales (sparse) para cada estado 'z'
    % spdiags crea una matriz sparse a partir de sus diagonales
    A1=spdiags(Y(:,1),0,I,I)+spdiags(X(2:I,1),-1,I,I)+spdiags([0;Z(1:I-1,1)],1,I,I);
    A2=spdiags(Y(:,2),0,I,I)+spdiags(X(2:I,2),-1,I,I)+spdiags([0;Z(1:I-1,2)],1,I,I);

    % Combina los bloques (movimiento en 'a') con la matriz de salto (movimiento en 'z')
    A = [A1,sparse(I,I);sparse(I,I),A2] + Aswitch; % Matriz completa (2I x 2I) 
    
    % Comprobación: Las filas de 'A' deben sumar 0 
    if max(abs(sum(A,2)))>10^(-9)
       disp('Matriz de Transición Impropia')
       break
    end
    
    % --- 5.4. RESOLVER EL SISTEMA (MÉTODO IMPLÍCITO) ---
    % La ecuación HJB discretizada en forma implícita es:
    % (v^{n+1} - v^n) / Delta + rho * v^{n+1} = u^n + A^n * v^{n+1} 
    % Reordenando:
    % [(1/Delta + rho)I - A^n] * v^{n+1} = u^n + v^n / Delta 
    % Llamamos B = [(1/Delta + rho)I - A^n] y b = u^n + v^n / Delta
    
    B = (1/Delta + rho)*speye(2*I) - A; % Matriz B
    
    % Apila (stack) los vectores para resolver el sistema lineal
    u_stacked = [u(:,1);u(:,2)];
    V_stacked = [V(:,1);V(:,2)]; % v^n
    
    b = u_stacked + V_stacked/Delta; % Vector b
    
    % Resuelve el sistema lineal B * v^{n+1} = b
    % Esta es la operación clave de la iteración implícita 
    V_stacked = B\b; % v^{n+1}
    
    % Des-apila (unstack) el vector de resultado para la próxima iteración
    V = [V_stacked(1:I),V_stacked(I+1:2*I)];
    
    % --- 5.5. Comprobación de Convergencia ---
    Vchange = V - v; % Diferencia entre v^{n+1} y v^n
    v = V; % Actualiza la función de valor para la próxima iteración

    dist(n) = max(max(abs(Vchange)));
    if dist(n)<crit
        disp('Función de Valor Convergida, Iteración = ')
        disp(n)
        break
    end
end
toc; % Fin del cronómetro para la HJB

% -------------------------------------------------------------------------
% 6. ECUACIÓN DE FOKKER-PLANCK (KOLMOGOROV FORWARD)
% -------------------------------------------------------------------------
% Ahora que 'A' ha convergido, encontramos la distribución estacionaria 'g'
% resolviendo A^T * g = 0 

AT = A'; % El operador de KF es la transpuesta (adjunta) del operador de HJB 
b = zeros(2*I,1); % El lado derecho es cero

% El sistema A^T * g = 0 es singular (las columnas de A^T suman 0).
% Necesitamos imponer la condición de que la distribución suma 1.
% El apéndice describe un "arreglo sucio" ("dirty fix"):
% reemplazar una fila de A^T por (0,...0,1,0,...0) y fijar b(i_fix) a un valor.
i_fix = 1;
b(i_fix)=.1; % Un valor arbitrario, será normalizado después
row = [zeros(1,i_fix-1),1,zeros(1,2*I-i_fix)];
AT(i_fix,:) = row; % Reemplaza la fila i_fix de A^T

% Resuelve el sistema lineal modificado A^T * g = b
gg = AT\b;

% Normaliza la distribución 'gg' para que la masa total sea 1
% La integral es la suma de (g_i * da) 
g_sum = gg'*ones(2*I,1)*da;
gg = gg./g_sum;

% Des-apila (unstack) el vector de distribución
g = [gg(1:I),gg(I+1:2*I)];

% Comprueba que la masa de cada tipo de ingreso suma a sus probabilidades estacionarias
check1 = g(:,1)'*ones(I,1)*da; % Masa total tipo 1 (debería ser la_2 / (la_1+la_2))
check2 = g(:,2)'*ones(I,1)*da; % Masa total tipo 2 (debería ser la_1 / (la_1+la_2))

% -------------------------------------------------------------------------
% 7. CÁLCULO DE AGREGADOS Y GRÁFICOS
% -------------------------------------------------------------------------

% Vuelve a calcular la política de ahorro (drift) final
adot = zz + r.*aa - c;
set(gca,'FontSize',14)
plot(a,adot,a,zeros(1,I),'--','LineWidth',2)
grid
xlabel('a')
ylabel('s_i(a)')
xlim([amin amax])

% Calcula la oferta de activos agregados S(r) 
S = g(:,1)'*a*da + g(:,2)'*a*da;


% --- Aproximación cerca de la restricción de endeudamiento ---
% Esto comprueba analíticamente la Proposición 1 del paper principal 
% s_1(a) \approx -\sqrt(2*\nu_1) * \sqrt(a - \underline{a}) 

% u'(c_1(\underline{a}))
u1 = (z1+r*amin)^(-s); 
% u'(c_2(\underline{a}))
u2 = c(1,2)^(-s);
% u''(c_1(\underline{a}))
u11 = -s*(z1+r*amin)^(-s-1);

% Fórmula para \nu_1 (nu) 
% \nu_1 = [ (rho-r)u'(c1) + la1(u'(c1) - u'(c2)) ] / (-u''(c1))
nu = sqrt(-2*((rho-r)*u1 + la1*(u1 - u2))/u11);

% Función de ahorro aproximada
s_approx = -nu*(a-amin).^(1/2);

% --- Gráficos ---
amax1 = 1; % Límite para el gráfico
set(gca,'FontSize',14)
h1 = plot(a,adot,a,s_approx,'-.',a,zeros(1,I),'--','LineWidth',2)
legend(h1,'s_1(a)','s_2(a)','Approximate s_1(a)','Location','SouthWest')
grid
xlabel('a')
ylabel('s_i(a)')
xlim([amin amax1])
% print -depsc HJB_stateconstraint.eps % (Comentado para no guardar archivo)

figure % Nueva figura para las densidades
set(gca,'FontSize',14)
h1 = plot(a,g,'LineWidth',2)
legend(h1,'g_1(a)','g_2(a)')
grid
xlabel('a')
ylabel('g_i(a)')
xlim([amin amax1])
% print -depsc densities.eps % (Comentado para no guardar archivo)