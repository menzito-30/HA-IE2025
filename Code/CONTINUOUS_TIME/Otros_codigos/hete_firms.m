% --- Inicialización ---
% Gracias a Peter Tian por sugerir este problema
clear all; close all; clc;
tic; % Inicia un cronómetro para medir el tiempo de ejecución

% --- Parámetros del Modelo ---
rho = 0.05; % Tasa de descuento (parámetro rho en la HJB)

% --- Definición del Proceso Estocástico (Ornstein-Uhlenbeck para log(z)) ---
% El código modela dlog(z) = -the*log(z)dt + sig*dW
% La distribución estacionaria es log(z) ~ N(0,Var) donde Var = sig2 / (2*the)
% 'the' es 'theta' (parámetro de reversión a la media), 'sig2' es la varianza del término de ruido
Var = 0.026^2;  % Varianza del proceso log(z)
zmean = exp(Var/2); % Media de z (dado que log(z) es normal)
Corr = 0.859; % Autocorrelación de log(z)
the = -log(Corr); % Parámetro de reversión a la media (theta)
sig2 = 2*the*Var; % Varianza del ruido (sigma^2)

% --- Grilla para z (Productividad) ---
J=40; % Número de puntos en la grilla de z
zmin = zmean*0.6; % Valor mínimo de z
zmax = zmean*1.4; % Valor máximo de z
z = linspace(zmin,zmax,J); % Vector de grilla para z
dz = (zmax-zmin)/(J-1); % Tamaño del paso en la grilla de z
dz2 = dz^2; % Tamaño del paso al cuadrado

% --- Grilla para k (Capital) ---
I=100; % Número de puntos en la grilla de k
kmin = 1; % Valor mínimo de k
kmax = 100; % Valor máximo de k
k = linspace(kmin,kmax,I)'; % Vector de grilla para k (como columna)
dk = (kmax-kmin)/(I-1); % Tamaño del paso en la grilla de k

% --- Creación de Matrices de Estado ---
% Se crean matrices 2D (I x J) para k y z que cubren todo el espacio de estados
kk = k*ones(1,J); % Matriz de capital (cada columna es 'k')
zz = ones(I,1)*z; % Matriz de productividad (cada fila es 'z')

% --- Parámetros del Proceso de Difusión para z ---
% Se aplica el Lema de Ito a z = exp(log(z)) para obtener dz = mu(z)dt + sigma(z)dW
% Esto coincide con la Ecuación (28) en HACT_Numerical_Appendix.pdf
mu = (-the*log(z) + sig2/2).*z; % Deriva (drift) del proceso z, mu(z)
s2 = sig2*(z.^2); % Varianza del proceso z, sigma^2(z)

% --- Parámetros de Iteración (Método Implícito) ---
maxit= 20; % Máximo de iteraciones
crit = 10^(-6); % Criterio de convergencia
Delta = 1000; % Tamaño del paso de tiempo para la iteración (Delta grande = método implícito)

% --- Parámetros Económicos ---
alpha= 0.5;     % Curvatura de la función de producción (k^alpha)
theta=2.7;    % Parámetro de coste de ajuste cuadrático. NOTA: el coste es 0.5*theta*(x/k)^2*k
delta=0.025;  % Tasa de depreciación
F = zz.*kk.^alpha; % Función de producción (matriz I x J)

% --- Pre-asignación de Matrices ---
Vkf = zeros(I,J); % Derivada de V respecto a k (forward)
Vkb = zeros(I,J); % Derivada de V respecto a k (backward)
Vzf = zeros(I,J); % Derivada de V respecto a z (forward)
Vzb = zeros(I,J); % Derivada de V respecto a z (backward)
Vzz = zeros(I,J); % Segunda derivada de V respecto a z
c = zeros(I,J);  % Variable de control (inversión 'x' en el paper)

% --- CONSTRUCCIÓN DE LA MATRIZ Bswitch (Transición para z) ---
% Esta matriz aproxima los términos de la HJB relacionados con la difusión de z:
% mu(z)*V_z + 0.5*sigma(z)^2*V_zz
% Se usan diferencias finitas (upwind para V_z, central para V_zz)
%
% chi: coeficiente para v_{i,j-1}
chi =  - min(mu,0)/dz + s2/(2*dz2);
% yy: coeficiente para v_{i,j}
yy =  min(mu,0)/dz - max(mu,0)/dz - s2/dz2;
% zeta: coeficiente para v_{i,j+1}
zeta = max(mu,0)/dz + s2/(2*dz2);

% 'Bswitch' es la parte de la matriz A total que maneja los saltos en z.
% Es la matriz 'C' en la Ecuación (36) y Figura 6 del apéndice

% Diagonal superior (para transiciones j -> j+1)
updiag=zeros(I,1); % Relleno inicial
for j=1:J
    updiag=[updiag;repmat(zeta(j),I,1)];
end

% Diagonal central (para transiciones j -> j)
% Se ajustan los bordes para reflejar
centdiag=repmat(chi(1)+yy(1),I,1);
for j=2:J-1
    centdiag=[centdiag;repmat(yy(j),I,1)];
end
centdiag=[centdiag;repmat(yy(J)+zeta(J),I,1)];

% Diagonal inferior (para transiciones j -> j-1)
lowdiag=repmat(chi(2),I,1);
for j=3:J
    lowdiag=[lowdiag;repmat(chi(j),I,1)];
end

% Construye la matriz 'Bswitch' (sparse)
% Las diagonales están separadas por I, ya que el vector de estado está apilado
% (k_1...k_I para z_1, luego k_1...k_I para z_2, etc.)
Bswitch=spdiags(centdiag,0,I*J,I*J)+spdiags(lowdiag,-I,I*J,I*J)+spdiags(updiag,I,I*J,I*J);

% --- Iteración de la Función de Valor (VFI) ---

% GUESS INICIAL
% El valor de estar en estado estacionario (x = delta*k) para siempre
% Flujo = F - delta*k - 0.5*theta*(delta*k/k)^2*k
v0 =(F - delta*kk - 0.5*theta*delta^2*kk)/rho; %
v = v0; % Asigna el guess inicial

plot(k,v) % Grafica el guess inicial

for n=1:maxit % Inicia el bucle de VFI
    V = v; % Almacena la función de valor de la iteración anterior
    
    % --- Esquema Upwind para k ---
    
    % Derivada hacia adelante (forward difference)
    Vkf(1:I-1,:) = (V(2:I,:)-V(1:I-1,:))/dk;
    % Condición de borde en k_max: se impone V_k = 1 + theta*delta
    % Esto implica x/k = delta (estado estacionario)
    Vkf(I,:) =   (1+theta*delta)*ones(1,J); 
    
    % Derivada hacia atrás (backward difference)
    Vkb(2:I,:) = (V(2:I,:)-V(1:I-1,:))/dk;
    % Condición de borde en k_min: se impone V_k = 1 + theta*delta
    Vkb(1,:) =   (1+theta*delta)*ones(1,J); 
    
    % --- Decisiones de Inversión (x) ---
    % Se deriva de la FOC: V_k = 1 + theta*(x/k)
    % (Nota: esto es por el 0.5*theta en la función de coste)
    
    % Inversión ('x') con derivada forward
    xf = (Vkf-1)/theta.*kk; % x = k * (V_k - 1) / theta
    sf = xf-delta.*kk; % Ahorro neto (drift de k): x - delta*k
    Hf = F - xf - 0.5*theta*(xf./kk).^2.*kk + Vkf.*sf; % Hamiltoniano
    
    % Inversión ('x') con derivada backward
    xb = (Vkb-1)/theta.*kk;
    sb =  xb-delta.*kk; 
    Hb = F - xb - 0.5*theta*(xb./kk).^2.*kk + Vkb.*sb;
    
    % Inversión en estado estacionario (para el caso I0)
    x0 =delta*kk;
    
    % Lógica del Upwind Scheme
    % Determina qué derivada usar basado en el signo del drift (s)
    Ineither = (1-(sf>0)) .* (1-(sb<0)); % Drift está entre 0 (s_f <= 0 <= s_b)
    Iunique = (sb<0).*(1-(sf>0)) + (1-(sb<0)).*(sf>0); % Solo una dirección es válida
    Iboth = (sb<0).*(sf>0); % Ambas son válidas (no debería pasar si V es cóncava)
    Ib = Iunique.*(sb<0) + Iboth.*(Hb>=Hf); % Casos para usar backward
    If = Iunique.*(sf>0) + Iboth.*(Hf>=Hb); % Casos para usar forward
    I0 = Ineither; % Caso para usar s=0 (x=delta*k)

    % Selecciona la inversión óptima 'x'
    x = xf.*If + xb.*Ib + x0.*I0;
    % Calcula el flujo de beneficios (payoff)
    profits = F - x - 0.5*theta*(x./kk).^2.*kk;
    
    % --- CONSTRUCCIÓN DE LA MATRIZ A (Transición para k) ---
    % Esta es la matriz A_tilde del apéndice
    
    % Coeficiente para v_{i-1,j}
    X = - sb.*Ib/dk;
    % Coeficiente para v_{i,j}
    Y = - sf.*If/dk + sb.*Ib/dk;
    % Coeficiente para v_{i+1,j}
    Z = sf.*If/dk;
    
    % Construcción de diagonales para la matriz AA (transiciones en k)
    updiag=0; 
    for j=1:J
        updiag=[updiag;Z(1:I-1,j);0]; % Diagonal superior
    end
    
    centdiag=reshape(Y,I*J,1); % Diagonal central
    
    lowdiag=X(2:I,1);
    for j=2:J
        lowdiag=[lowdiag;0;X(2:I,j)]; % Diagonal inferior
    end
    
    % Matriz AA (sparse) para transiciones en k
    AA=spdiags(centdiag,0,I*J,I*J)+spdiags([updiag;0],1,I*J,I*J)+spdiags([lowdiag;0],-1,I*J,I*J);
    
    % Matriz A total (transiciones en k y z)
    A = AA + Bswitch;
    
    % Chequeo: las filas de A deben sumar (casi) cero
    if max(abs(sum(A,2)))>10^(-12)
        disp('Error: La matriz de transición no suma cero.')
        break
    end
    
    % --- Resolución del Sistema Implícito ---
    % Se resuelve B * V^{n+1} = b
    % donde B = (1/Delta + rho)I - A
    % y b = profits + V^{n} / Delta
    
    B = (1/Delta + rho)*speye(I*J) - A;
    
    profits_stacked = reshape(profits,I*J,1); % Apila vector de beneficios
    V_stacked = reshape(V,I*J,1); % Apila función de valor
    
    b = profits_stacked + V_stacked/Delta; % Construye el vector b

    V_stacked = B\b; % Resuelve el sistema lineal para V^{n+1}
    
    V = reshape(V_stacked,I,J); % Des-apila la nueva V
    
    % --- Chequeo de Convergencia ---
    Vchange = V - v; % Diferencia con la iteración anterior
    v = V; % Actualiza la función de valor

    dist(n) = max(max(abs(Vchange))); % Almacena la máxima diferencia
    if dist(n)<crit
        disp('Función de Valor Convergió, Iteración = ')
        disp(n)
        break % Termina el bucle si converge
    end
end
toc; % Detiene el cronómetro

% --- Gráficos Finales ---
figure;
plot(k,v) % Grafica la función de valor convergida
title('Función de Valor Convergida')

figure;
plot(dist) % Grafica la distancia (error) en cada iteración
title('Convergencia')

kdot = x - delta.*kk; % Calcula la política de ahorro neto (k_dot)
figure;
plot(k,kdot,k,zeros(I,1),'--') % Grafica la política de ahorro
title('Política de Ahorro (kdot) vs k')
legend('kdot (z=z_{min})', 'kdot (z=z_{max})', 'kdot=0')