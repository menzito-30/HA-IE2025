clear all; clc; close all;

% -------------------------------------------------------------------------
% 1. DEFINICIÓN DE PARÁMETROS DEL MODELO
% -------------------------------------------------------------------------

s = 2;            % Coeficiente de aversión relativa al riesgo (\gamma) 
rho = 0.05;       % Tasa de descuento subjetiva
r = 0.03;         % Tasa de interés
z1 = .1;          % Productividad baja
z2 = .2;          % Productividad alta
z = [z1,z2];      % Vector de estados de productividad
la1 = 1.5;        % Intensidad de Poisson para pasar de z1 a z2 (\lambda_1)
la2 = 1;          % Intensidad de Poisson para pasar de z2 a z1 (\lambda_2)
la = [la1,la2];   % Vector de intensidades

frisch = 1/2;     % Elasticidad de Frisch de la oferta laboral (\varphi)
w = 1;            % Salario (fijo)

% -------------------------------------------------------------------------
% 2. CREACIÓN DE LA GRILLA (REJILLA)
% -------------------------------------------------------------------------
% Discretización del espacio de estados 'a'

I = 500;             % Número de puntos en la grilla de activos 'a'
amin = -0.15;      % Límite inferior de activos (\underline{a})
amax = 3;          % Límite superior de activos
a = linspace(amin,amax,I)'; % Vector columna de la grilla de activos
da = (amax-amin)/(I-1);  % Distancia entre puntos de la grilla (\Delta a)

% Matrices auxiliares para cálculos vectorizados
aa = [a,a];        % Repite la grilla de activos para ambos estados de 'z'
zz = ones(I,1)*z;  % Matriz (I x 2) con los valores de 'z'

% -------------------------------------------------------------------------
% 3. PARÁMETROS DE ITERACIÓN NUMÉRICA
% -------------------------------------------------------------------------

maxit= 20;         % Número máximo de iteraciones para la función de valor
crit = 10^(-6);    % Criterio de convergencia (tolerancia)
Delta = 1000;      % Tamaño del paso para el método implícito

% Inicialización de matrices
dVf = zeros(I,2); % Derivada aproximada hacia adelante (Forward)
dVb = zeros(I,2); % Derivada aproximada hacia atrás (Backward)
c = zeros(I,2);   % Política de consumo

% Matriz (sparse) que captura los saltos entre estados de productividad (z1 <-> z2)
Aswitch = [-speye(I)*la(1),speye(I)*la(1);speye(I)*la(2),-speye(I)*la(2)];

% -------------------------------------------------------------------------
% 4. CÁLCULO PREVIO: OFERTA LABORAL CON AHORRO CERO (l0)
% -------------------------------------------------------------------------
%  Lo más "complicado" es manejar los
% puntos donde el ahorro (drift) es cero. En esos puntos, se deben
% cumplir dos condiciones simultáneamente:
% 1. Ahorro cero: c = w*z*l + r*a 
% 2. FOC Intratemporal: l^(1/phi) * c^gamma = w*z 
%
% Este bucle resuelve ese sistema no lineal (usando 'fzero') para CADA
% punto 'a' de la grilla ANTES de empezar la iteración de la HJB.
% 'l0' almacenará la oferta laboral que genera ahorro cero en cada (a, z).

% Opciones para el solver 'fzero' (para no mostrar la salida)
options=optimset('Display','off');
% Conjetura inicial para el solver 'fzero'
x0 = (w*z1)^(frisch*(1-s)/(1+s*frisch));
   
tic; % Inicia cronómetro para este cálculo previo
for i=1:I
   % Resuelve para el estado z1 en el punto a(i)
   params = [a(i),z1,w,r,s,frisch]; % Parámetros para la función 'lab_solve'
   myfun = @(l) lab_solve(l,params); % 'lab_solve' debe implementar la ec. (36)
   [l01,fval,exitflag] = fzero(myfun,x0,options); % Encuentra 'l'

   % Resuelve para el estado z2 en el punto a(i)
   params = [a(i),z2,w,r,s,frisch];
   myfun = @(l) lab_solve(l,params);
   [l02,fval,exitflag] = fzero(myfun,x0,options); % Encuentra 'l'
   
   % Almacena la oferta laboral de ahorro cero
   l0(i,:)=[l01,l02];
end
toc % Fin del cronómetro

% -------------------------------------------------------------------------
% 5. CONJETURA INICIAL (INITIAL GUESS) PARA LA FUNCIÓN DE VALOR
% -------------------------------------------------------------------------
% (Nota: esta conjetura usa el 'l0' de la frontera para todos los 'a',
% es solo un punto de partida).
v0(:,1) = (w*z(1).*l0(1,1) + r.*a).^(1-s)/(1-s)/rho;
v0(:,2) = (w*z(2).*l0(1,2) + r.*a).^(1-s)/(1-s)/rho;

% Almacena los valores de 'l0' en las fronteras
lmin = l0(1,:); % l0 en amin (se usará para la condición de contorno)
lmax = l0(I,:); % l0 en amax

v = v0; % 'v' es la función de valor que se actualizará en cada iteración

% -------------------------------------------------------------------------
% 6. BUCLE PRINCIPAL: ITERACIÓN DE LA FUNCIÓN DE VALOR (HJB)
% -------------------------------------------------------------------------
%maxit = 1;
for n=1:maxit
    V = v; % Guarda la función de valor de la iteración anterior (v^n)
    V_n(:,:,n)=V;
    
    % --- 6.1. Aproximación de derivadas (Diferencias Finitas) ---
    
    % Derivada hacia adelante: v'(a_i) \approx (v_{i+1} - v_i) / da
    dVf(1:I-1,:) = (V(2:I,:)-V(1:I-1,:))/da;
    % Condición de contorno en amax (usando l0 en amax)
    dVf(I,:) = (w*z.*lmax + r.*amax).^(-s); 
    
    % Derivada hacia atrás: v'(a_i) \approx (v_i - v_{i-1}) / da
    dVb(2:I,:) = (V(2:I,:)-V(1:I-1,:))/da;
    % Condición de restricción estatal en amin 
    % v'(\underline{a}) = u_c(\underline{c}, \underline{l})
    % Se impone con igualdad usando c = w*z*l + r*a 
    dVb(1,:) = (w*z.*lmin + r.*amin).^(-s); 
       
    % --- 6.2. ESQUEMA "UPWIND" ---
    
    % 1. Decisiones con derivada hacia adelante (forward, v'_F)
    % De FOC u_c = v' => c = (v')^(-1/s)
    cf = dVf.^(-1/s);
    % De FOC intratemporal => l^(1/phi) = w*z / c^gamma
    % Pero como u_c = c^(-gamma) = v', esto es l = (v' * w * z)^phi 
    lf = (dVf.*w.*zz).^frisch;
    % Ahorro (drift): s_j(a) = w*z_j*l_j + r*a - c_j 
    ssf = w*zz.*lf + r.*aa - cf;
    
    % 2. Decisiones con derivada hacia atrás (backward, v'_B)
    cb = dVb.^(-1/s);
    lb = (dVb.*w.*zz).^frisch;
    ssb = w*zz.*lb + r.*aa - cb;
    
    % 3. Decisiones con ahorro cero (drift = 0)
    % Usamos el 'l0' pre-calculado
    c0 = w*zz.*l0 + r.*aa; % Consumo consistente con ahorro cero 
    dV0 = c0.^(-s);        % v' consistente con c0
    
    % Indicadores para elegir la derivada correcta (esquema upwind)
    Ib = ssb < 0; % Usar 'backward' si el drift (ahorro) es negativo
    If = (ssf > 0).*(1-Ib); % Usar 'forward' si el drift es positivo
    I0 = (1-If-Ib); % Usar 'cero' (l0, c0) si el drift cambia de signo
    
    % Selecciona las políticas óptimas de consumo y trabajo
    c = cf.*If + cb.*Ib + c0.*I0;
    l = lf.*If + lb.*Ib + l0.*I0;
    
    % Calcula la utilidad con la desutilidad del trabajo 
    u = c.^(1-s)/(1-s) - l.^(1+1/frisch)/(1+1/frisch);
     
    % --- 6.3. CONSTRUCCIÓN DE LA MATRIZ DE TRANSICIÓN 'A' ---
    
    X = -Ib.*ssb/da; % Coeficiente para v_{i-1} (viene de 'backward')
    Y = -If.*ssf/da + Ib.*ssb/da; % Coeficiente para v_i
    Z = If.*ssf/da; % Coeficiente para v_{i+1} (viene de 'forward')
    
    A1=spdiags(Y(:,1),0,I,I)+spdiags(X(2:I,1),-1,I,I)+spdiags([0;Z(1:I-1,1)],1,I,I);
    A2=spdiags(Y(:,2),0,I,I)+spdiags(X(2:I,2),-1,I,I)+spdiags([0;Z(1:I-1,2)],1,I,I);

    A = [A1,sparse(I,I);sparse(I,I),A2] + Aswitch;
    
    % Comprobación: Las filas de 'A' deben sumar 0
    if max(abs(sum(A,2)))>10^(-9)
       disp('Matriz de Transición Impropia')
       break
    end
    
    % --- 6.4. RESOLVER EL SISTEMA (MÉTODO IMPLÍCITO) ---
    % Reordenamos la HJB:
    % [(1/Delta + rho)I - A^n] * v^{n+1} = u^n + v^n / Delta
    % B * v^{n+1} = b
    
    B = (1/Delta + rho)*speye(2*I) - A; % Matriz B
    
    % Apila (stack) los vectores para resolver el sistema lineal
    u_stacked = [u(:,1);u(:,2)];
    V_stacked = [V(:,1);V(:,2)]; % v^n
    
    b = u_stacked + V_stacked/Delta; % Vector b
    
    V_stacked = B\b; % Resuelve para v^{n+1}
      
    % Des-apila (unstack) el vector de resultado
    V = [V_stacked(1:I),V_stacked(I+1:2*I)];
    
    % --- 6.5. Comprobación de Convergencia ---
    Vchange = V - v;
    v = V; % Actualiza la función de valor

    dist(n) = max(max(abs(Vchange)));
    if dist(n)<crit
        disp('Función de Valor Convergida, Iteración = ')
        disp(n)
        break
    end
end
toc; % Fin del cronómetro para la HJB

%%%%%%%%%%%%%%%%%%%%%%%%%%
% FOKKER-PLANCK EQUATION %
%%%%%%%%%%%%%%%%%%%%%%%%%%
AT = A';
b = zeros(2*I,1);

%need to fix one value, otherwise matrix is singular
i_fix = 1;
b(i_fix)=.1;
row = [zeros(1,i_fix-1),1,zeros(1,2*I-i_fix)];
AT(i_fix,:) = row;

%Solve linear system
gg = AT\b;
g_sum = gg'*ones(2*I,1)*da;
gg = gg./g_sum;

g = [gg(1:I),gg(I+1:2*I)];

check1 = g(:,1)'*ones(I,1)*da;
check2 = g(:,2)'*ones(I,1)*da;



% -------------------------------------------------------------------------
% 7. GRÁFICOS DE RESULTADOS
% -------------------------------------------------------------------------

% Política de ahorro final
adot = w.*zz.*l + r.*aa - c;

% Gráfico de la política de ahorro (s_j(a))
figure(1)
plot(a,adot,a,zeros(1,I),'--','LineWidth',2)
set(gca,'FontSize',14)
legend('s_1(a)','s_2(a)','Location','SouthWest')
grid
xlabel('Wealth, a')
ylabel('Saving, s_i(a)')
% print -depsc labor_s.eps % (Comentado para no guardar archivo)

% Gráfico de la función de valor (v_j(a))
figure(2)
plot(a,v,'LineWidth',2)
set(gca,'FontSize',14)
legend('v_1(a)','v_2(a)','Location','SouthWest')
grid
xlabel('Wealth, a')
ylabel('Value Function, v_i(a)')
% print -depsc labor_v.eps % (Comentado para no guardar archivo)

% Gráfico de la política de oferta laboral (l_j(a))
figure(3)
plot(a,l,'LineWidth',2)
set(gca,'FontSize',14)
legend('l_1(a)','l_2(a)','Location','SouthWest')
grid
xlabel('Wealth, a')
ylabel('Labor Supply, l_i(a)')
% print -depsc labor_l.eps % (Comentado para no guardar archivo)

% Gráfico de la política de consumo (c_j(a))
figure(4)
plot(a,c,'LineWidth',2)
set(gca,'FontSize',14)
legend('c_1(a)','c_2(a)','Location','SouthWest')
grid
xlabel('Wealth, a')
ylabel('Consumption, c_i(a)')
% print -depsc labor_c.eps % (Comentado para no guardar archivo)


figure(5) % Nueva figura
plot(a,g,'LineWidth',2)
set(gca,'FontSize',14)
legend('g_1(a)','g_2(a)','Location','NorthEast')
grid
xlabel('Wealth, a')
ylabel('Density, g_i(a)')