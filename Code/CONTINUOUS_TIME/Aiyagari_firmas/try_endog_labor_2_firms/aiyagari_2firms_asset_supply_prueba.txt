%% AIYAGARI MODEL WITH ENDOGENOUS LABOR SUPPLY - TWO SECTORS (FORMAL & INFORMAL)
% Based on Achdou, Han, Lasry, Lions, and Moll (2017)
% Extension with two sectors: Formal and Informal
% demora aprox 6 minutos en correr!

clear all; clc; close all;

tic; % Start timer

% =========================================================================
% 1. PARAMETERS
% =========================================================================

% --- Household Parameters ---
ga = 2;           % Relative risk aversion (γ)
rho = 0.05;       % Subjective discount rate (ρ)
Frisch = 0.5;     % Frisch elasticity of labor supply (φ)

% --- Productivity States (z) ---
z1 = 0.2;         % Low productivity
z2 = 2*z1;        % High productivity
z = [z1, z2];     % Productivity vector
la1 = 1;          % Poisson intensity z1 -> z2
la2 = 1;          % Poisson intensity z2 -> z1
la = [la1, la2];

% Average productivity (for initial labor supply normalization)
z_ave = (z1*la2 + z2*la1)/(la1 + la2);

% --- Firm Parameters ---
Aprod = 0.3;      % Total factor productivity (A)
al = 1/3;         % Capital share (α)
d = 0.05;         % Depreciation rate (δ)

% --- Sector-Specific Parameters ---
AF = Aprod;       % TFP Formal sector
AI = 0.8*Aprod;   % TFP Informal sector (lower productivity)

% =========================================================================
% 2. GRIDS
% =========================================================================

% --- Asset Grid ---
I = 500;          % Number of grid points
amin = 0;         % Borrowing constraint
amax = 20;        % Maximum assets
a = linspace(amin, amax, I)'; % Asset grid
da = (amax - amin)/(I - 1);   % Grid spacing

% Auxiliary matrices
aa = [a, a];      % Replicate asset grid for both z states
zz = ones(I,1)*z; % Matrix (I x 2) with z values

% --- HJB Iteration Parameters ---
maxit = 100;      % Maximum iterations
crit = 10^(-6);   % Convergence criterion
Delta = 1000;     % Step size for implicit method

% Initialize matrices
dVf = zeros(I, 2); % Forward derivative
dVb = zeros(I, 2); % Backward derivative
c = zeros(I, 2);   % Consumption policy
lf = zeros(I, 2);  % Labor supply formal sector
li = zeros(I, 2);  % Labor supply informal sector

% Transition matrix for productivity shocks
Aswitch = [-speye(I)*la(1), speye(I)*la(1); 
           speye(I)*la(2), -speye(I)*la(2)];

% --- Price Grid ---
Ir = 40;          % Number of interest rate points
rmin = -0.045;    % Minimum interest rate
rmax = 0.045;     % Maximum interest rate
r_grid = linspace(rmin, rmax, Ir);

% =========================================================================
% 3. INITIAL GUESS (for first r)
% =========================================================================
r = r_grid(1);

% --- Initial Firm Prices (assuming single firm for initial guess) ---
KD_guess = (al*Aprod/(r + d))^(1/(1 - al))*z_ave;
wF = (1 - al)*Aprod*(KD_guess/z_ave)^al;
wI = 0.7*wF;  % Informal wage is lower (initial guess)

% --- Initial Value Function Guess ---
% Assume constant labor supply of 1 for initial guess
v0(:,1) = (wF*z(1) + max(r,0.01)*a).^(1 - ga)/(1 - ga)/rho;
v0(:,2) = (wF*z(2) + max(r,0.01)*a).^(1 - ga)/(1 - ga)/rho;

% =========================================================================
% 4. OUTER LOOP: ITERATE OVER PRICES
% =========================================================================

% Storage
S = zeros(Ir, 1);        % Asset supply
KF_demand = zeros(Ir, 1); % Capital demand formal
KI_demand = zeros(Ir, 1); % Capital demand informal
wF_r = zeros(Ir, 1);     % Wage formal vs r
wI_r = zeros(Ir, 1);     % Wage informal vs r
LF_s = zeros(Ir, 1);     % Labor supply formal
LI_s = zeros(Ir, 1);     % Labor supply informal

for ir = 1:Ir
    
    r = r_grid(ir);
    
    fprintf('\n=== Iteration %d/%d, r = %.4f ===\n', ir, Ir, r)
    
    % --- 4.1. Firm Prices (given r) ---
    % For now, assume equal capital allocation K = KF + KI
    % We'll iterate on wage ratio wI/wF
    
    % Initial guess: equal capital-labor ratios
    kappa = (al*AF/(r + d))^(1/(1-al));
    wF = (1 - al)*AF*kappa^al;
    wI = (1 - al)*AI*kappa^al * (AI/AF)^(1/(1-al));
    
    wF_r(ir) = wF;
    wI_r(ir) = wI;
    
    % --- 4.2. Warm Start ---
    if ir > 1
        v0 = V_r(:,:,ir - 1);
    end
    
    v = v0;
    
    % Precompute boundary labor for state constraint
    lf_min = zeros(1,2);
    li_min = zeros(1,2);
    c_min = zeros(1,2);
    dV_min = zeros(1,2);
    
    for j = 1:2
        income_min = @(l) wF*z(j)*l(1) + wI*z(j)*l(2) + r*amin;
        eq_min = @(l) [l(1)^(1/Frisch) * (income_min(l))^ga - wF*z(j);
                       l(2)^(1/Frisch) * (income_min(l))^ga - wI*z(j)];
        l_init = [0.5; 0.5];
        options = optimoptions('fsolve', 'Display', 'off', 'TolFun', 1e-10);
        l_sol = fsolve(eq_min, l_init, options);
        lf_min(j) = max(l_sol(1), 0);
        li_min(j) = max(l_sol(2), 0);
        c_min(j) = income_min(l_sol);
        dV_min(j) = c_min(j)^(-ga);
    end
    
    % =====================================================================
    % 5. INNER LOOP: SOLVE HJB (given r, wF, wI)
    % =====================================================================
    
    for n = 1:maxit
        V = v;
        
        % --- 5.1. Finite Difference Approximations ---
        % Forward derivative
        dVf(1:I-1,:) = (V(2:I,:) - V(1:I-1,:))/da;
        % Backward derivative
        dVb(2:I,:) = (V(2:I,:) - V(1:I-1,:))/da;
        
        % --- 5.2. Boundary Conditions ---
        % At amax: Use FOC with labor supply at boundary
        for j = 1:2
            params_upper = [amax, z(j), wF, wI, r, ga, Frisch];
            l_upper = fsolve(@(l) lab_solve_dual(l, params_upper), [1; 1], ...
                optimoptions('fsolve', 'Display', 'off', 'TolFun', 1e-10));
            income_upper = wF*z(j)*l_upper(1) + wI*z(j)*l_upper(2) + r*amax;
            dVf(I,j) = income_upper^(-ga);
        end
        
        % At amin: Enforce state constraint
        dVb(1,:) = max( (V(2,:) - V(1,:))/da, dV_min );
        
        % --- 5.3. SOLVE FOR LABOR SUPPLY AND CONSUMPTION ---
        % Forward direction
        dVf_pos = max(dVf, 1e-10);
        cf = real(dVf_pos .^ (-1/ga));
        lf_f = zeros(I, 2);
        li_f = zeros(I, 2);
        
        for j = 1:2
            for i = 1:I
                params = [a(i), z(j), wF, wI, r, ga, Frisch];
                l_guess = [0.5; 0.5];
                try
                    l_sol = fsolve(@(l) lab_solve_dual(l, params), l_guess, ...
                        optimoptions('fsolve', 'Display', 'off', 'TolFun', 1e-10));
                    lf_f(i,j) = max(l_sol(1), 0);
                    li_f(i,j) = max(l_sol(2), 0);
                catch
                    lf_f(i,j) = 0.5;
                    li_f(i,j) = 0.5;
                end
            end
        end
        ssf = wF*zz.*lf_f + wI*zz.*li_f + r*aa - cf;
        
        % Backward direction
        dVb_pos = max(dVb, 1e-10);
        cb = real(dVb_pos .^ (-1/ga));
        lf_b = zeros(I, 2);
        li_b = zeros(I, 2);
        
        for j = 1:2
            for i = 1:I
                params = [a(i), z(j), wF, wI, r, ga, Frisch];
                l_guess = [0.5; 0.5];
                try
                    l_sol = fsolve(@(l) lab_solve_dual(l, params), l_guess, ...
                        optimoptions('fsolve', 'Display', 'off', 'TolFun', 1e-10));
                    lf_b(i,j) = max(l_sol(1), 0);
                    li_b(i,j) = max(l_sol(2), 0);
                catch
                    lf_b(i,j) = 0.5;
                    li_b(i,j) = 0.5;
                end
            end
        end
        ssb = wF*zz.*lf_b + wI*zz.*li_b + r*aa - cb;
        
        % Zero drift
        c0 = zeros(I, 2);
        lf_0 = zeros(I, 2);
        li_0 = zeros(I, 2);
        
        for j = 1:2
            for i = 1:I
                params = [a(i), z(j), wF, wI, r, ga, Frisch];
                l_guess = [0.5; 0.5];
                try
                    l_sol = fsolve(@(l) lab_solve_dual(l, params), l_guess, ...
                        optimoptions('fsolve', 'Display', 'off', 'TolFun', 1e-10));
                    lf_0(i,j) = max(l_sol(1), 0);
                    li_0(i,j) = max(l_sol(2), 0);
                catch
                    lf_0(i,j) = 0.5;
                    li_0(i,j) = 0.5;
                end
                income_0 = wF*z(j)*lf_0(i,j) + wI*z(j)*li_0(i,j) + r*a(i);
                c0(i,j) = income_0;
            end
        end
        dV0 = max(c0 .^ (-ga), 1e-10);
        
        % Upwind indicators
        If = ssf > 0;
        Ib = ssb < 0 & ~If;
        I0 = ~(If | Ib);
        
        % Construct upwind derivative and policies
        dV_Upwind = dVf_pos.*If + dVb_pos.*Ib + dV0.*I0;
        c = cf.*If + cb.*Ib + c0.*I0;
        lf = lf_f.*If + lf_b.*Ib + lf_0.*I0;
        li = li_f.*If + li_b.*Ib + li_0.*I0;
        
        % Utility
        u = c.^(1 - ga)/(1 - ga) - lf.^(1 + 1/Frisch)/(1 + 1/Frisch) ...
            - li.^(1 + 1/Frisch)/(1 + 1/Frisch);
        
        % --- 5.4. CONSTRUCT TRANSITION MATRIX A ---
        ss_Upwind = ssf .* If + ssb .* Ib + 0 .* I0;
        X = -min(ss_Upwind, 0) / da;
        Y = -max(ss_Upwind, 0) / da + min(ss_Upwind, 0) / da;
        Z = max(ss_Upwind, 0) / da;
        
        % Construct A1 and A2
        A1 = spdiags(Y(:,1), 0, I, I) + spdiags([X(2:I,1); 0], -1, I, I) + ...
             spdiags([0; Z(1:I-1,1)], 1, I, I);
        A2 = spdiags(Y(:,2), 0, I, I) + spdiags([X(2:I,2); 0], -1, I, I) + ...
             spdiags([0; Z(1:I-1,2)], 1, I, I);
        A = [A1, sparse(I,I); sparse(I,I), A2] + Aswitch;
        
        % Verify row sums
        row_sum = sum(A, 2);
        if max(abs(row_sum)) > 1e-9
            fprintf('Warning: Improper Transition Matrix at iteration %d\n', n);
        end
        
        % --- 5.5. SOLVE IMPLICIT SYSTEM ---
        B = (1/Delta + rho)*speye(2*I) - A;
        
        u_stacked = [u(:,1); u(:,2)];
        V_stacked = [V(:,1); V(:,2)];
        
        b = u_stacked + V_stacked/Delta;
        V_stacked = B\b;
        
        V = [V_stacked(1:I), V_stacked(I+1:2*I)];
        
        % --- 5.6. Check Convergence ---
        Vchange = V - v;
        v = V;
        dist(n) = max(max(abs(Vchange)));
        
        if dist(n) < crit
            fprintf('Value Function Converged at iteration %d for r = %.4f\n', n, r)
            break
        end
    end
    
    % =====================================================================
    % 6. KOLMOGOROV FORWARD EQUATION
    % =====================================================================
    
    AT = A';
    b_KF = zeros(2*I, 1);
    
    % Fix first element
    i_fix = 1;
    b_KF(i_fix) = 0.1;
    row = [zeros(1, i_fix-1), 1, zeros(1, 2*I - i_fix)];
    AT(i_fix,:) = row;
    
    gg = AT\b_KF;
    g_sum = gg'*ones(2*I,1)*da;
    gg = gg/g_sum;
    
    g = [gg(1:I), gg(I+1:2*I)];
    
    % --- 6.1. Store Results ---
    g_r(:,:,ir) = g;
    adot(:,:,ir) = wF*zz.*lf + wI*zz.*li + r*aa - c;
    V_r(:,:,ir) = V;
    dV_r(:,:,ir) = dV_Upwind;
    c_r(:,:,ir) = c;
    lf_r(:,:,ir) = lf;
    li_r(:,:,ir) = li;
    
    % --- KEY CALCULATIONS ---
    % Asset supply
    S(ir) = g(:,1)'*a*da + g(:,2)'*a*da;
    
    % Labor supply by sector
    LF_s(ir) = da * (g(:,1)' * (z(1) * lf(:,1)) + g(:,2)' * (z(2) * lf(:,2)));
    LI_s(ir) = da * (g(:,1)' * (z(1) * li(:,1)) + g(:,2)' * (z(2) * li(:,2)));
    
    % Capital demand (from firm FOCs)
    % For simplicity, assume equal capital-labor ratios initially
    % In full equilibrium, need to iterate on capital allocation
    total_L = LF_s(ir) + LI_s(ir);
    KF_demand(ir) = (al*AF/(r + d))^(1/(1-al)) * LF_s(ir);
    KI_demand(ir) = (al*AI/(r + d))^(1/(1-al)) * LI_s(ir);
    
    % Display progress
    if mod(ir, 10) == 0
        fprintf('r = %.4f: S = %.4f, KF = %.4f, KI = %.4f, LF = %.4f, LI = %.4f\n', ...
            r, S(ir), KF_demand(ir), KI_demand(ir), LF_s(ir), LI_s(ir))
    end
    
end

toc; % End timer

% =========================================================================
% 7. FIND EQUILIBRIUM
% =========================================================================

K_demand = KF_demand + KI_demand;
excess_supply = S - K_demand;

% Find equilibrium interest rate
[~, idx_eq] = min(abs(excess_supply));
r_eq = r_grid(idx_eq);
K_eq = S(idx_eq);
wF_eq = wF_r(idx_eq);
wI_eq = wI_r(idx_eq);
LF_eq = LF_s(idx_eq);
LI_eq = LI_s(idx_eq);

fprintf('\n=== EQUILIBRIUM ===\n')
fprintf('r* = %.4f\n', r_eq)
fprintf('K* = %.4f\n', K_eq)
fprintf('wF* = %.4f\n', wF_eq)
fprintf('wI* = %.4f\n', wI_eq)
fprintf('LF* = %.4f\n', LF_eq)
fprintf('LI* = %.4f\n', LI_eq)
fprintf('Wage gap (wF/wI) = %.2f\n', wF_eq/wI_eq)

% =========================================================================
% 8. PLOT RESULTS
% =========================================================================

figure('Position', [100, 100, 1200, 800])

% Asset supply and demand
subplot(2,3,1)
plot(S, r_grid, 'b-', 'LineWidth', 2, 'DisplayName', 'S(r) - Asset Supply')
hold on
plot(K_demand, r_grid, 'r-', 'LineWidth', 2, 'DisplayName', 'K^D(r) - Capital Demand')
plot(K_eq, r_eq, 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'g', 'DisplayName', 'Equilibrium')
xlabel('Capital/Assets', 'FontSize', 12)
ylabel('Interest Rate r', 'FontSize', 12)
title('Asset Market Equilibrium', 'FontSize', 14)
legend('show', 'Location', 'best')
grid on

% Labor supply by sector
subplot(2,3,2)
plot(r_grid, LF_s, 'b-', 'LineWidth', 2, 'DisplayName', 'L^F(r) - Formal')
hold on
plot(r_grid, LI_s, 'r-', 'LineWidth', 2, 'DisplayName', 'L^I(r) - Informal')
plot(r_eq, LF_eq, 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b')
plot(r_eq, LI_eq, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r')
xlabel('Interest Rate r', 'FontSize', 12)
ylabel('Labor Supply', 'FontSize', 12)
title('Labor Supply by Sector', 'FontSize', 14)
legend('show', 'Location', 'best')
grid on

% Wage gap
subplot(2,3,3)
plot(r_grid, wF_r./wI_r, 'k-', 'LineWidth', 2)
hold on
plot(r_eq, wF_eq/wI_eq, 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'g')
xlabel('Interest Rate r', 'FontSize', 12)
ylabel('Wage Ratio w^F/w^I', 'FontSize', 12)
title('Formal-Informal Wage Gap', 'FontSize', 14)
grid on

% Distribution
subplot(2,3,4)
g_eq = g_r(:,:,idx_eq);
plot(a, g_eq(:,1), 'b-', 'LineWidth', 2, 'DisplayName', 'z = z_1')
hold on
plot(a, g_eq(:,2), 'r-', 'LineWidth', 2, 'DisplayName', 'z = z_2')
xlabel('Assets a', 'FontSize', 12)
ylabel('Density g(a,z)', 'FontSize', 12)
title('Wealth Distribution', 'FontSize', 14)
legend('show')
grid on

% Labor policies at equilibrium
subplot(2,3,5)
lf_eq = lf_r(:,:,idx_eq);
li_eq = li_r(:,:,idx_eq);
plot(a, lf_eq(:,1), 'b-', 'LineWidth', 2, 'DisplayName', 'l^F, z=z_1')
hold on
plot(a, lf_eq(:,2), 'b--', 'LineWidth', 2, 'DisplayName', 'l^F, z=z_2')
plot(a, li_eq(:,1), 'r-', 'LineWidth', 2, 'DisplayName', 'l^I, z=z_1')
plot(a, li_eq(:,2), 'r--', 'LineWidth', 2, 'DisplayName', 'l^I, z=z_2')
xlabel('Assets a', 'FontSize', 12)
ylabel('Labor Supply', 'FontSize', 12)
title('Labor Supply Policies', 'FontSize', 14)
legend('show', 'Location', 'best')
grid on

% Consumption policy
subplot(2,3,6)
c_eq = c_r(:,:,idx_eq);
plot(a, c_eq(:,1), 'b-', 'LineWidth', 2, 'DisplayName', 'z = z_1')
hold on
plot(a, c_eq(:,2), 'r-', 'LineWidth', 2, 'DisplayName', 'z = z_2')
xlabel('Assets a', 'FontSize', 12)
ylabel('Consumption c(a,z)', 'FontSize', 12)
title('Consumption Policy', 'FontSize', 14)
legend('show')
grid on

% Save figure
print('-dpng', 'aiyagari_two_sectors_results.png', '-r300')

% =========================================================================
% 9. SAVE RESULTS
% =========================================================================

save('aiyagari_two_sectors_results.mat', ...
     'r_grid', 'r_eq', 'K_eq', 'wF_eq', 'wI_eq', 'LF_eq', 'LI_eq', ...
     'S', 'K_demand', 'KF_demand', 'KI_demand', 'LF_s', 'LI_s', ...
     'V_r', 'g_r', 'c_r', 'lf_r', 'li_r', 'adot', ...
     'a', 'z', 'ga', 'Frisch', 'rho', 'al', 'd', 'AF', 'AI')

fprintf('\n=== SIMULATION COMPLETE ===\n')
fprintf('Results saved to: aiyagari_two_sectors_results.mat\n')
fprintf('Figure saved to: aiyagari_two_sectors_results.png\n')