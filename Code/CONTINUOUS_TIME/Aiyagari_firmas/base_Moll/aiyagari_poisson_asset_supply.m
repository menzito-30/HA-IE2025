%Optimized for speed by SeHyoun Ahn

clear all; clc; close all;

tic; % Inicia el cronómetro

% -------------------------------------------------------------------------
% 1. PARÁMETROS DEL MODELO (HOGARES Y FIRMAS)
% -------------------------------------------------------------------------

% --- Parámetros de los Hogares ---
s = 2;            % Coeficiente de aversión relativa al riesgo (\gamma)
rho = 0.05;       % Tasa de descuento subjetiva
z1 = .2;          % Productividad baja
z2 = 2*z1;        % Productividad alta
z = [z1,z2];      % Vector de estados de productividad
la1 = 1;          % Tasa de Poisson z1 -> z2 (\lambda_1)
la2 = 1;          % Tasa de Poisson z2 -> z1 (\lambda_2)
la = [la1,la2];
% Trabajo agregado (L), se asume igual a la media estacionaria de 'z'
z_ave = (z1*la2+z2*la1)/(la1+la2);

% --- Parámetros de las Firmas (Producción) ---
% Se asume una función de producción Cobb-Douglas: Y = A * K^al * L^(1-al)
Aprod = 0.3;      % Productividad total de los factores (A)
al = 1/3;         % Participación del capital en la producción (\alpha)
d = 0.05;         % Tasa de depreciación del capital (\delta)


% -------------------------------------------------------------------------
% 2. CREACIÓN DE LAS GRILLAS (REJILLAS)
% -------------------------------------------------------------------------

% --- Grilla de Activos (Capital 'a') ---
I=500;
amin = 0;         % Límite inferior (Restricción de no-endeudamiento, \underline{a}=0)
amax = 20;
a = linspace(amin,amax,I)'; % Vector columna de la grilla de activos
da = (amax-amin)/(I-1);  % Distancia entre puntos de la grilla (\Delta a)

% Matrices auxiliares para cálculos vectorizados
aa = [a,a];        % Repite la grilla de activos para ambos estados de 'z'
zz = ones(I,1)*z;  % Matriz (I x 2) con los valores de 'z'

% --- Parámetros de Iteración (HJB) ---
maxit= 100;        % Número máximo de iteraciones para la HJB
crit = 10^(-6);    % Criterio de convergencia (tolerancia) para la HJB
Delta = 1000;      % Tamaño del paso para el método implícito

% Inicialización de matrices
dVf = zeros(I,2);  % Derivada HJB hacia adelante
dVb = zeros(I,2);  % Derivada HJB hacia atrás
c = zeros(I,2);    % Política de consumo

% Matriz (sparse) que captura los saltos entre estados de productividad (z1 <-> z2)
Aswitch = [-speye(I)*la(1),speye(I)*la(1);speye(I)*la(2),-speye(I)*la(2)];

% --- Grilla de Tasa de Interés ---
% El objetivo de este script es trazar la curva S(r)
Ir = 100;         % Número de puntos en la grilla de 'r'
rmin = -0.0499;   % Tasa de interés mínima (cerca de -delta)
rmax = 0.049;     % Tasa de interés máxima (cerca de rho)
r_grid = linspace(rmin,rmax,Ir); % Grilla de tasas de interés a evaluar

% -------------------------------------------------------------------------
% 3. CONJETURA INICIAL (PARA LA PRIMERA 'r' DE LA GRILLA)
% -------------------------------------------------------------------------
r = r_grid(1); % Empieza con la tasa de interés más baja

% --- Problema de la Firma (para la 'r' inicial) ---
% Dado 'r', la firma demanda capital (KD) y trabajo (L=z_ave)
% FOC Capital: r = F_K - d = al*Aprod*(KD/z_ave)^(al-1) - d
KD = (al*Aprod/(r+d))^(1/(1-al))*z_ave; % Demanda de capital para la r inicial
% FOC Trabajo: w = F_L = (1-al)*Aprod*(KD/z_ave)^al
w = (1-al)*Aprod*(KD/z_ave)^al; % Salario para la r inicial

% --- Conjetura Inicial HJB (para la 'r' inicial) ---
% Se usa una fórmula basada en el consumo perpetuo del ingreso
% (Se usa max(r, 0.01) para evitar problemas si r es muy negativo)
v0(:,1) = (w*z(1) + max(r,0.01).*a).^(1-s)/(1-s)/rho;
v0(:,2) = (w*z(2) + max(r,0.01).*a).^(1-s)/(1-s)/rho;

% -------------------------------------------------------------------------
% 4. BUCLE EXTERNO: ITERAR SOBRE LA GRILLA DE TASAS DE INTERÉS 'r'
% -------------------------------------------------------------------------
% Para cada 'r' en 'r_grid', se resuelve el problema del hogar (HJB)
% y se calcula la oferta agregada de capital (S(r)) usando la KF.

for ir=1:Ir;

r = r_grid(ir); % Fija la tasa de interés para esta iteración

% --- 4.1. Resolver Precios de la Firma (Dado 'r') ---
% Se calcula la Demanda de Capital (KD) y el Salario (w) que son
% consistentes con la tasa de interés 'r' actual.
KD(ir) = (al*Aprod/(r+d))^(1/(1-al))*z_ave; % Demanda de Capital K_D(r)
w = (1-al)*Aprod*(KD(ir)/z_ave)^al; % Salario w(r)
w_r(ir)=w; % Almacena el salario para esta 'r'

% --- 4.2. Conjetura Inicial HJB (Warm Start) ---
% Si no es la primera iteración (ir > 1), usa la función de valor
% convergida de la 'r' anterior como conjetura inicial. Acelera mucho.
if ir>1
v0 = V_r(:,:,ir-1);
end

v = v0; % 'v' es la función de valor actual que se actualizará

% -------------------------------------------------------------------------
% 5. BUCLE INTERNO: RESOLVER HJB (Dado 'r' y 'w')
% -------------------------------------------------------------------------
% Este bucle encuentra la función de valor v(a,z) y las políticas c(a,z), s(a,z)
% para la tasa de interés 'r' y el salario 'w' actuales.
for n=1:maxit
    V = v; % Guarda el valor de la iteración anterior (v^n)
    V_n(:,:,n)=V; % Almacena historial (opcional)
    
    % --- 5.1. Aproximación de derivadas (Diferencias Finitas) ---
    % Derivada hacia adelante
    dVf(1:I-1,:) = (V(2:I,:)-V(1:I-1,:))/da;
    dVf(I,:) = (w*z + r.*amax).^(-s); % Condición de contorno en amax
    % Derivada hacia atrás
    dVb(2:I,:) = (V(2:I,:)-V(1:I-1,:))/da;
    dVb(1,:) = (w*z + r.*amin).^(-s); % Condición de contorno en amin (a=0)
    
    I_concave = dVb > dVf; % Indicador de concavidad (para depuración)
    
    % --- 5.2. ESQUEMA "UPWIND" ---
    % 1. Ahorro (drift) con derivada hacia adelante (forward)
    cf = dVf.^(-1/s); % c = (v')^{-1/s}
    ssf = w*zz + r.*aa - cf; % s_j(a) = w*z_j + r*a - c_j(a)
    % 2. Ahorro (drift) con derivada hacia atrás (backward)
    cb = dVb.^(-1/s);
    ssb = w*zz + r.*aa - cb;
    % 3. Ahorro nulo (drift = 0)
    c0 = w*zz + r.*aa; % c = w*z + r*a
    dV0 = c0.^(-s);    % v' correspondiente
    
    % Indicadores para elegir la derivada correcta (esquema upwind)
    If = ssf > 0; % Usar 'forward' si el drift es positivo
    Ib = ssb < 0; % Usar 'backward' si el drift es negativo
    I0 = (1-If-Ib); % Usar 'cero' (c0, dV0) si el drift cambia de signo
    
    % Construye la derivada upwind y la política de consumo final
    dV_Upwind = dVf.*If + dVb.*Ib + dV0.*I0;
    c = dV_Upwind.^(-1/s);
    % Calcula la utilidad
    u = c.^(1-s)/(1-s);
    
    % --- 5.3. CONSTRUCCIÓN DE LA MATRIZ DE TRANSICIÓN 'A' ---
    % 'A' es el generador infinitesimal discretizado de la HJB
    X = - min(ssb,0)/da; % Coeficiente para v_{i-1,j}
    Y = - max(ssf,0)/da + min(ssb,0)/da; % Coeficiente para v_{i,j}
    Z = max(ssf,0)/da; % Coeficiente para v_{i+1,j}
    
    A1=spdiags(Y(:,1),0,I,I)+spdiags(X(2:I,1),-1,I,I)+spdiags([0;Z(1:I-1,1)],1,I,I);
    A2=spdiags(Y(:,2),0,I,I)+spdiags(X(2:I,2),-1,I,I)+spdiags([0;Z(1:I-1,2)],1,I,I);
    % Matriz 'A' completa (incluye saltos de 'z')
    A = [A1,sparse(I,I);sparse(I,I),A2] + Aswitch;
    
    % Comprobación: Las filas de 'A' deben sumar 0
    if max(abs(sum(A,2)))>10^(-9)
       disp('Improper Transition Matrix')
       break
    end
    
    % --- 5.4. RESOLVER EL SISTEMA IMPLÍCITO ---
    % La HJB discretizada es: B * v^{n+1} = b
    B = (1/Delta + rho)*speye(2*I) - A; % Matriz B
    
    % Apila (stack) los vectores
    u_stacked = [u(:,1);u(:,2)];
    V_stacked = [V(:,1);V(:,2)]; % v^n
    
    b = u_stacked + V_stacked/Delta; % Vector b
    
    V_stacked = B\b; % Resuelve para v^{n+1}
    
    % Des-apila (unstack) el vector de resultado
    V = [V_stacked(1:I),V_stacked(I+1:2*I)];
    
    % --- 5.5. Comprobación de Convergencia HJB ---
    Vchange = V - v;
    v = V; % Actualiza la función de valor
    dist(n) = max(max(abs(Vchange)));
    if dist(n)<crit
        disp('Value Function Converged, Iteration = ')
        disp(n)
        break
    end
end % Fin del bucle interno (HJB)
toc; % Fin del cronómetro para ESTA iteración de HJB

% -------------------------------------------------------------------------
% 6. FOKKER-PLANCK (KF) (Dado 'A' convergido para este 'r')
% -------------------------------------------------------------------------
% Ahora se resuelve A^T * g = 0 para encontrar la distribución estacionaria
% asociada a las políticas óptimas para la 'r' actual.

AT = A'; % El operador de KF es la transpuesta (adjunta) del operador de HJB
b = zeros(2*I,1); % El lado derecho es cero

% "Dirty fix" para la singularidad: fija g(1)=0.1 y reemplaza la 1ra ecuación
i_fix = 1;
b(i_fix)=.1;
row = [zeros(1,i_fix-1),1,zeros(1,2*I-i_fix)];
AT(i_fix,:) = row;

% Resuelve el sistema lineal modificado A^T * g = b
gg = AT\b;
% Normaliza la distribución 'gg' para que la masa total (integral) sea 1
g_sum = gg'*ones(2*I,1)*da; % Integral = suma de (g_i * da)
gg = gg./g_sum;

% Des-apila (unstack) el vector de distribución
g = [gg(1:I),gg(I+1:2*I)];

% (Opcional) Comprueba que las masas por tipo suman a las prob. estacionarias
check1 = g(:,1)'*ones(I,1)*da;
check2 = g(:,2)'*ones(I,1)*da;

% --- 6.1. Almacenar resultados para esta iteración de 'r' ---
g_r(:,:,ir) = g; % Almacena la distribución estacionaria para r(ir)
adot(:,:,ir) = w*zz + r.*aa - c; % Almacena la política de ahorro para r(ir)
V_r(:,:,ir) = V; % Almacena la función de valor para r(ir)
dV_r(:,:,ir) = dV_Upwind; % Almacena la derivada v' para r(ir)
c_r(:,:,ir) = c; % Almacena la política de consumo para r(ir)

% --- ¡CÁLCULO CLAVE PARA ESTE SCRIPT! ---
% Calcular la Oferta Agregada de Activos (Capital) S(r) para esta 'r'
% K_S = \int a * g(a,z | r) da dz
S(ir) = g(:,1)'*a*da + g(:,2)'*a*da; % Almacena S(r(ir))

end % Fin del bucle externo (sobre 'r')

% -------------------------------------------------------------------------
% 7. GRÁFICO DE EQUILIBRIO DE MERCADO (OFERTA Y DEMANDA)
% -------------------------------------------------------------------------
% Prepara variables auxiliares para el gráfico
Smax = max(S); % Límite superior para el eje K
amin1 = amin-0.02; % Límite inferior para el eje K
aaa = linspace(amin1,Smax,Ir); % Rango para el eje K
rrr = linspace(-0.06,0.06,Ir); % Rango para el eje r
% Recalcula la curva de Demanda de Capital (K_D) para el rango del gráfico
KD = (al*Aprod./(max(rrr+d,0))).^(1/(1-al))*z_ave;

% Crea el gráfico
set(gca,'FontSize',14)
% Grafica S(r) vs K_D(r)
plot(S,r_grid,KD,rrr,zeros(1,Ir)+amin,rrr,'--',aaa,ones(1,Ir)*rho,'--',aaa,ones(1,Ir)*(-d),'--','LineWidth',2)
% Añade etiquetas y líneas guía
text(0.05,0.052,'$r = \rho$','FontSize',16,'interpreter','latex')
text(0.05,-0.054,'$r = -\delta$','FontSize',16,'interpreter','latex')
text(0.1,0.02,'$S(r)$','FontSize',16,'interpreter','latex') % Curva de Oferta de Capital (Hogares)
text(0.29,0.02,'$r=F_K(K,1)-\delta$','FontSize',16,'interpreter','latex') % Curva de Demanda de Capital (Firmas)
text(0.01,0.035,'$a=\underline{a}$','FontSize',16,'interpreter','latex') % Línea de restricción a=0
ylabel('$r$','FontSize',16,'interpreter','latex') % Eje Y: Tasa de interés
xlabel('$K$','FontSize',16,'interpreter','latex') % Eje X: Capital
% Ajusta los límites de los ejes
ylim([-0.06 0.06])
xlim([amin1 0.6])
% Guarda el gráfico (comentado)
% print -depsc aiyagari_asset_supply.eps

% El punto donde S(r) y K_D(r) se cruzan es el equilibrio general estacionario.