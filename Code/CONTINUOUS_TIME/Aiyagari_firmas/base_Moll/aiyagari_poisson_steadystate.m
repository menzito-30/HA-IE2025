%Optimized for speed by SeHyoun Ahn

clear all; clc; close all;

tic; % Inicia el cronómetro

% -------------------------------------------------------------------------
% 1. DEFINICIÓN DE PARÁMETROS DEL MODELO
% -------------------------------------------------------------------------

% --- Parámetros de los Hogares ---
ga = 2;            % Coeficiente de aversión relativa al riesgo (\gamma)
rho = 0.05;        % Tasa de descuento subjetiva
d = 0.05;          % Tasa de depreciación del capital (\delta)

% --- Parámetros de las Firmas (Producción) ---
al = 1/3;          % Participación del capital (\alpha) en Cobb-Douglas Y = A*K^al*L^(1-al)
Aprod = 0.1;       % Productividad total de los factores (A)

% --- Parámetros del Proceso Estocástico (Ingreso) ---
z1 = 1;            % Productividad baja
z2 = 2*z1;         % Productividad alta
z = [z1,z2];       % Vector de estados de productividad
la1 = 1/3;         % Tasa de Poisson z1 -> z2 (\lambda_1)
la2 = 1/3;         % Tasa de Poisson z2 -> z1 (\lambda_2)
la = [la1,la2];
% Trabajo agregado (L), se asume igual a la media estacionaria de 'z'
z_ave = (z1*la2 + z2*la1)/(la1 + la2);

% -------------------------------------------------------------------------
% 2. CREACIÓN DE LA GRILLA (REJILLA) DE ACTIVOS
% -------------------------------------------------------------------------

I= 1000;           % Número de puntos en la grilla de activos 'a'
amin = 0;          % Límite inferior (Restricción de no-endeudamiento, \underline{a}=0)
amax = 20;         % Límite superior de activos
% Comentario del código original: Sugerencia para reescalar 'amax' si 'Aprod' cambia
% amax = 20*(Aprod/0.1)^(1/(1-al));

a = linspace(amin,amax,I)'; % Vector columna de la grilla de activos
da = (amax-amin)/(I-1);  % Distancia entre puntos de la grilla (\Delta a)

% Matrices auxiliares para cálculos vectorizados
aa = [a,a];        % Repite la grilla de activos para ambos estados de 'z'
zz = ones(I,1)*z;  % Matriz (I x 2) con los valores de 'z'

% -------------------------------------------------------------------------
% 3. PARÁMETROS DE ITERACIÓN NUMÉRICA
% -------------------------------------------------------------------------

% --- Parámetros HJB (Bucle Interno) ---
maxit= 100;        % Número máximo de iteraciones para la HJB
crit = 10^(-6);    % Criterio de convergencia (tolerancia) para la HJB
Delta = 1000;      % Tamaño del paso para el método implícito

% Inicialización de matrices
dVf = zeros(I,2);  % Derivada HJB hacia adelante
dVb = zeros(I,2);  % Derivada HJB hacia atrás
c = zeros(I,2);    % Política de consumo

% Matriz (sparse) que captura los saltos entre estados de productividad (z1 <-> z2)
Aswitch = [-speye(I)*la(1),speye(I)*la(1);speye(I)*la(2),-speye(I)*la(2)];

% --- Parámetros Equilibrio (Bucle Externo) ---
Ir = 40;           % Máximas iteraciones para encontrar 'r' (Equilibrio General)
crit_S = 10^(-5);  % Criterio de convergencia para el mercado de capitales (S=0)

% (Valores iniciales para 'r' y 'w' - serán sobreescritos)
rmax = 0.049;
r = 0.04;
w = 0.05;

% --- Parámetros de Bisección (Búsqueda de 'r') ---
r0 = 0.03;         % (No se usa, pero está definido)
rmin = 0.01;       % Límite inferior para la búsqueda de 'r'
rmax = 0.99*rho;   % Límite superior para la búsqueda de 'r' (r < rho)

% -------------------------------------------------------------------------
% 4. BUCLE EXTERNO: BÚSQUEDA DEL EQUILIBRIO GENERAL (BISECCIÓN)
% -------------------------------------------------------------------------
% Este bucle busca la tasa de interés 'r' que equilibra el mercado de capitales,
% es decir, donde la Oferta de Capital (KS) = Demanda de Capital (KD).

for ir=1:Ir;

% Almacena el historial de la búsqueda
r_r(ir)=r;
rmin_r(ir)=rmin;
rmax_r(ir)=rmax;

% --- 4.1. Resolver Precios de la Firma (Dado 'r') ---
% A partir de las FOC de la firma:
% FOC Capital: r = F_K - d = al*Aprod*(L/K)^(1-al) - d
% Despejando K (que es KD, la Demanda de Capital):
KD(ir) = (al*Aprod/(r + d))^(1/(1-al))*z_ave;
% FOC Trabajo: w = F_L = (1-al)*Aprod*(K/L)^al
w = (1-al)*Aprod*KD(ir).^al*z_ave^(-al); % Salario de equilibrio

% Chequeo de que la restricción de endeudamiento es válida (natural)
if w*z(1) + r*amin < 0
    disp('CAREFUL: borrowing constraint too loose')
end

% --- 4.2. Conjetura Inicial HJB (Warm Start) ---
% Para la primera iteración (ir=1), se usa una fórmula analítica
% Para ir>1, se usa la función de valor (V) convergida de la iteración anterior
if ir==1
    v0(:,1) = (w*z(1) + r.*a).^(1-ga)/(1-ga)/rho;
    v0(:,2) = (w*z(2) + r.*a).^(1-ga)/(1-ga)/rho;
else
    v0 = V_r(:,:,ir-1); % "Warm start"
end

v = v0; % 'v' es la función de valor actual que se actualizará

% -------------------------------------------------------------------------
% 5. BUCLE INTERNO: RESOLVER HJB (Dado 'r' y 'w')
% -------------------------------------------------------------------------
for n=1:maxit
    V = v; % Guarda el valor de la iteración anterior (v^n)
    V_n(:,:,n)=V; % Almacena historial (opcional)
    
    % --- 5.1. Aproximación de derivadas (Diferencias Finitas) ---
    % Derivada hacia adelante: v'(a_i) \approx (v_{i+1} - v_i) / da
    dVf(1:I-1,:) = (V(2:I,:)-V(1:I-1,:))/da;
    dVf(I,:) = (w*z + r.*amax).^(-ga); % Condición de contorno en amax
    
    % Derivada hacia atrás: v'(a_i) \approx (v_i - v_{i-1}) / da
    dVb(2:I,:) = (V(2:I,:)-V(1:I-1,:))/da;
    dVb(1,:) = (w*z + r.*amin).^(-ga); % Condición de contorno en amin
    
    % --- 5.2. ESQUEMA "UPWIND" ---
    % 1. Ahorro (drift) con derivada hacia adelante (forward)
    cf = dVf.^(-1/ga); % c = (u')^{-1}(v') = (v')^{-1/ga}
    ssf = w*zz + r.*aa - cf; % s_j(a) = w*z_j + r*a - c_j(a)
    
    % 2. Ahorro (drift) con derivada hacia atrás (backward)
    cb = dVb.^(-1/ga);
    ssb = w*zz + r.*aa - cb;
    
    % 3. Ahorro nulo (drift = 0)
    c0 = w*zz + r.*aa; % c = w*z + r*a
    
    % Indicadores para elegir la derivada correcta (esquema upwind)
    If = ssf > 0; % Usar 'forward' si el drift (ahorro) es positivo
    Ib = ssb < 0; % Usar 'backward' si el drift (ahorro) es negativo
    I0 = (1-If-Ib); % Usar 'cero' (c0) si el drift cambia de signo
    
    % Construye la política de consumo final
    c = cf.*If + cb.*Ib + c0.*I0;
    % Calcula la utilidad
    u = c.^(1-ga)/(1-ga);
    
    % --- 5.3. CONSTRUCCIÓN DE LA MATRIZ DE TRANSICIÓN 'A' ---
    % 'A' es el generador infinitesimal discretizado
    X = -min(ssb,0)/da; % Coeficiente para v_{i-1,j}
    Y = -max(ssf,0)/da + min(ssb,0)/da; % Coeficiente para v_{i,j}
    Z = max(ssf,0)/da; % Coeficiente para v_{i+1,j}
    
    A1=spdiags(Y(:,1),0,I,I)+spdiags(X(2:I,1),-1,I,I)+spdiags([0;Z(1:I-1,1)],1,I,I);
    A2=spdiags(Y(:,2),0,I,I)+spdiags(X(2:I,2),-1,I,I)+spdiags([0;Z(1:I-1,2)],1,I,I);

    % Matriz 'A' completa (incluye saltos de 'z')
    A = [A1,sparse(I,I);sparse(I,I),A2] + Aswitch;
    
    % Comprobación: Las filas de 'A' (el generador) deben sumar 0
    if max(abs(sum(A,2)))>10^(-9)
       disp('Improper Transition Matrix')
       %break
    end
    
    % --- 5.4. RESOLVER EL SISTEMA IMPLÍCITO ---
    % La HJB discretizada es: (1/Delta + rho)v^{n+1} - A v^{n+1} = u^n + v^n/Delta
    %                     o: B * v^{n+1} = b
    
    B = (1/Delta + rho)*speye(2*I) - A; % Matriz B
    
    % Apila (stack) los vectores para resolver el sistema lineal
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

AT = A'; % El operador de KF es la transpuesta (adjunta) del operador de HJB
b = zeros(2*I,1); % El lado derecho es cero

% "Dirty fix": El sistema A^T * g = 0 es singular (g=0 es solución).
% Se impone una condición (g_1 = 0.1) para encontrar una solución no trivial
% que luego se normaliza para que la masa total sea 1.
i_fix = 1;
b(i_fix)=.1;
row = [zeros(1,i_fix-1),1,zeros(1,2*I-i_fix)];
AT(i_fix,:) = row; % Reemplaza la primera fila de A^T

% Resuelve el sistema lineal modificado A^T * g = b
gg = AT\b;
% Normaliza la distribución 'gg' para que la masa total (integral) sea 1
g_sum = gg'*ones(2*I,1)*da; % Integral = suma de (g_i * da)
gg = gg./g_sum;

% Des-apila (unstack) el vector de distribución
g = [gg(1:I),gg(I+1:2*I)];

% (Opcional) Comprueba que la masa de cada tipo suma a sus probabilidades estacionarias
check1 = g(:,1)'*ones(I,1)*da;
check2 = g(:,2)'*ones(I,1)*da;

% --- 6.1. Almacenar resultados para esta iteración de 'r' ---
g_r(:,:,ir) = g; % Almacena la distribución
adot(:,:,ir) = w*zz + r.*aa - c; % Almacena la política de ahorro
V_r(:,:,ir) = V; % Almacena la función de valor

% -------------------------------------------------------------------------
% 7. ACTUALIZACIÓN DEL BUCLE DE EQUILIBRIO (BISECCIÓN)
% -------------------------------------------------------------------------

% Calcular la Oferta de Capital de los hogares K_S(r)
% K_S = \int a * g(a,z) da dz
KS(ir) = g(:,1)'*a*da + g(:,2)'*a*da; 
% Calcular el Exceso de Oferta de Capital: S(r) = K_S(r) - K_D(r)
S(ir) = KS(ir) - KD(ir); % (KD se calculó en el paso 4.1)

% --- Lógica de Bisección ---
if S(ir)>crit_S % Exceso de Oferta de capital (hogares ahorran demasiado)
    disp('Excess Supply')
    rmax = r; % El 'r' de equilibrio debe ser más bajo (para desincentivar el ahorro)
    r = 0.5*(r+rmin); % Prueba un 'r' a medio camino entre r y rmin
elseif S(ir)<-crit_S; % Exceso de Demanda de capital (hogares ahorran muy poco)
    disp('Excess Demand')
    rmin = r; % El 'r' de equilibrio debe ser más alto (para incentivar el ahorro)
    r = 0.5*(r+rmax); % Prueba un 'r' a medio camino entre r y rmax
elseif abs(S(ir))<crit_S; % ¡Equilibrio encontrado! (Exceso de Oferta/Demanda es ~0)
    display('Equilibrium Found, Interest rate =')
    disp(r)
    break % Termina el bucle externo 'ir'
end

end % Fin del bucle externo (sobre 'r')

% -------------------------------------------------------------------------
% 8. GRÁFICOS DE RESULTADOS DEL EQUILIBRIO
% -------------------------------------------------------------------------
% (El código restante es para crear y formatear los gráficos finales
% usando los resultados de la última iteración 'ir')

amax1 = 5;
amin1 = amin-0.1;

figure(1)
h1 = plot(a,adot(:,1,ir),'b',a,adot(:,2,ir),'r',linspace(amin1,amax1,I),zeros(1,I),'k--','LineWidth',2)
legend(h1,'s_1(a)','s_2(a)','Location','NorthEast')
text(-0.155,-.105,'$\underline{a}$','FontSize',16,'interpreter','latex')
line([amin amin], [-.1 .08],'Color','Black','LineStyle','--')
xlabel('Wealth, $a$','interpreter','latex')
ylabel('Savings, $s_i(a)$','interpreter','latex')
xlim([amin1 amax1])
ylim([-0.03 0.05])
set(gca,'FontSize',16)

figure(2)
h1 = plot(a,g_r(:,1,ir),'b',a,g_r(:,2,ir),'r','LineWidth',2)
legend(h1,'g_1(a)','g_2(a)')
text(-0.155,-.12,'$\underline{a}$','FontSize',16,'interpreter','latex')
line([amin amin], [0 max(max(g_r(:,:,ir)))],'Color','Black','LineStyle','--')
xlabel('Wealth, $a$','interpreter','latex')
ylabel('Densities, $g_i(a)$','interpreter','latex')
xlim([amin1 amax1])
%ylim([0 0.5])
set(gca,'FontSize',16)

