%% AIYAGARI MODEL WITH ENDOGENOUS LABOR SUPPLY - ASSET SUPPLY CURVE
% Based on Achdou, Han, Lasry, Lions, and Moll (2017)
% Extension with endogenous labor supply following Moll's labor_supply.pdf
% Author: Adapted for endogenous labor
% Date: 2025

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
ell = zeros(I, 2); % Labor supply policy

% Transition matrix for productivity shocks
Aswitch = [-speye(I)*la(1), speye(I)*la(1); 
           speye(I)*la(2), -speye(I)*la(2)];

% --- Interest Rate Grid ---
Ir = 100;         % Number of interest rate points
rmin = -0.0499;   % Minimum interest rate
rmax = 0.049;     % Maximum interest rate
r_grid = linspace(rmin, rmax, Ir);

% =========================================================================
% 3. INITIAL GUESS (for first r)
% =========================================================================
r = r_grid(1);

% --- Firm Problem (given r) ---
KD = (al*Aprod/(r + d))^(1/(1 - al))*z_ave;
w = (1 - al)*Aprod*(KD/z_ave)^al;

% --- Initial Value Function Guess ---
% Assume constant labor supply of 1 for initial guess
v0(:,1) = (w*z(1) + max(r,0.01)*a).^(1 - ga)/(1 - ga)/rho;
v0(:,2) = (w*z(2) + max(r,0.01)*a).^(1 - ga)/(1 - ga)/rho;

% =========================================================================
% 4. OUTER LOOP: ITERATE OVER INTEREST RATES
% =========================================================================

for ir = 1:Ir
    
    r = r_grid(ir);
    
    % --- 4.1. Firm Prices (given r) ---
    KD(ir) = (al*Aprod/(r + d))^(1/(1 - al))*z_ave;
    w = (1 - al)*Aprod*(KD(ir)/z_ave)^al;
    w_r(ir) = w;
    
    % --- 4.2. Warm Start ---
    if ir > 1
        v0 = V_r(:,:,ir - 1);
    end
    
    v = v0;
    
    % Precompute boundary labor for state constraint
    ell_min = zeros(1,2);
    c_min = zeros(1,2);
    dV_min = zeros(1,2);
    for j = 1:2
        income_min = @(l) w*z(j)*l + r*amin;
        eq_min = @(l) l^(1/Frisch) * (income_min(l))^ga - w*z(j);
        ell_min(j) = fzero(eq_min, 1, optimset('Display','off','TolFun',1e-10));
        c_min(j) = income_min(ell_min(j));
        dV_min(j) = c_min(j)^(-ga);
    end
    
    % =====================================================================
    % 5. INNER LOOP: SOLVE HJB (given r and w)
    % =====================================================================
    
    for n = 1:maxit
        V = v;
        % V_n(:,:,n) = V; % Optional: store history
        
        % --- 5.1. Finite Difference Approximations ---
        % Forward derivative
        dVf(1:I-1,:) = (V(2:I,:) - V(1:I-1,:))/da;
        % Backward derivative
        dVb(2:I,:) = (V(2:I,:) - V(1:I-1,:))/da;
        
        % --- 5.2. Boundary Conditions ---
        % At amax: Use FOC with labor supply at boundary
        dVf(I,:) = (w*z + r.*amax).^(-ga); % <-- Revertir a la BC estándar
        
        % At amin: Enforce state constraint
        dVb(1,:) = max( (V(2,:) - V(1,:))/da, dV_min );  % Ensure >= constraint
        
        % --- 5.3. UPWIND SCHEME WITH ENDOGENOUS LABOR (Revised) ---
        % Forward difference
        dVf_pos = max(dVf, 1e-10);  % Avoid negative derivatives
        cf = real(dVf_pos .^ (-1/ga));
        %ell_f = zeros(I, 2);
        for j = 1:2
            for i = 1:I
                ell_guess = 1;
                eq_f = @(l) l^(1/Frisch) * cf(i,j)^ga - w*z(j);
                ell_f(i,j) = fzero(eq_f, ell_guess, optimset('Display','off','TolFun',1e-10));
                if isnan(ell_f(i,j))  % Fallback if fzero fails
                    ell_f(i,j) = 1;
                end
            end
        end
        ssf = w*zz.*ell_f + r*aa - cf;
        
        % Backward difference
        dVb_pos = max(dVb, 1e-10);
        cb = real(dVb_pos .^ (-1/ga));
        %ell_b = zeros(I, 2);
        for j = 1:2
            for i = 1:I
                ell_guess = 1;
                eq_b = @(l) l^(1/Frisch) * cb(i,j)^ga - w*z(j);
                ell_b(i,j) = fzero(eq_b, ell_guess, optimset('Display','off','TolFun',1e-10));
                if isnan(ell_b(i,j))
                    ell_b(i,j) = 1;
                end
            end
        end
        ssb = w*zz.*ell_b + r*aa - cb;
        
        % Zero drift
        c0 = zeros(I, 2);
        ell_0 = zeros(I, 2);
        for j = 1:2
            for i = 1:I
                income_0 = @(l) w*z(j)*l + r*a(i);
                eq_0 = @(l) l^(1/Frisch) * (income_0(l))^ga - w*z(j);
                ell_0(i,j) = fzero(eq_0, 1, optimset('Display','off','TolFun',1e-10));
                if isnan(ell_0(i,j))
                    ell_0(i,j) = 1;
                end
                c0(i,j) = income_0(ell_0(i,j));
            end
        end
        dV0 = max(c0 .^ (-ga), 1e-10);  % Ensure positive
        
        % Upwind indicators
        If = ssf > 0;
        Ib = ssb < 0 & ~If;
        I0 = ~(If | Ib);
        
        % Construct upwind derivative and policies
        dV_Upwind = dVf_pos.*If + dVb_pos.*Ib + dV0.*I0;
        c = cf.*If + cb.*Ib + c0.*I0;
        ell = ell_f.*If + ell_b.*Ib + ell_0.*I0;
        
        % Utility
        u = c.^(1 - ga)/(1 - ga) - ell.^(1 + 1/Frisch)/(1 + 1/Frisch);
        
        % --- 5.4. CONSTRUCT TRANSITION MATRIX A (Revised) ---
        % Use upwind savings for consistency
        X = - min(ssb,0)/da; % Coeficiente para v_{i-1,j}
        Y = - max(ssf,0)/da + min(ssb,0)/da; % Coeficiente para v_{i,j}
        Z = max(ssf,0)/da; % Coeficiente para v_{i+1,j}
                
        % Construct A1 and A2
        A1 = spdiags(Y(:,1), 0, I, I) + spdiags([X(2:I,1); 0], -1, I, I) + ...
             spdiags([0; Z(1:I-1,1)], 1, I, I);
        A2 = spdiags(Y(:,2), 0, I, I) + spdiags([X(2:I,2); 0], -1, I, I) + ...
             spdiags([0; Z(1:I-1,2)], 1, I, I);
        A = [A1, sparse(I,I); sparse(I,I), A2] + Aswitch;
        
        % Verify and adjust row sums
        row_sum = sum(A, 2);
        if max(abs(row_sum)) > 1e-9
            disp('Improper Transition Matrix - Adjusting Y');
            Y(:,1) = Y(:,1) - row_sum(1:I) / I;
            Y(:,2) = Y(:,2) - row_sum(I+1:2*I) / I;
            % Reconstruct
            A1 = spdiags(Y(:,1), 0, I, I) + spdiags([X(2:I,1); 0], -1, I, I) + ...
                 spdiags([0; Z(1:I-1,1)], 1, I, I);
            A2 = spdiags(Y(:,2), 0, I, I) + spdiags([X(2:I,2); 0], -1, I, I) + ...
                 spdiags([0; Z(1:I-1,2)], 1, I, I);
            A = [A1, sparse(I,I); sparse(I,I), A2] + Aswitch;
            if max(abs(sum(A, 2))) > 1e-9
                disp('Adjustment failed - Check drift terms');
                break;
            end
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
    b = zeros(2*I, 1);
    
    % Fix first element
    i_fix = 1;
    b(i_fix) = 0.1;
    row = [zeros(1, i_fix-1), 1, zeros(1, 2*I - i_fix)];
    AT(i_fix,:) = row;
    
    gg = AT\b;
    g_sum = gg'*ones(2*I,1)*da;
    gg = gg/g_sum;
    
    g = [gg(1:I), gg(I+1:2*I)];
    
    % Checks
    check1 = g(:,1)'*ones(I,1)*da;
    check2 = g(:,2)'*ones(I,1)*da;
    
    % --- 6.1. Store Results ---
    g_r(:,:,ir) = g;
    adot(:,:,ir) = w*zz.*ell + r*aa - c;  % Savings policy
    V_r(:,:,ir) = V;
    dV_r(:,:,ir) = dV_Upwind;
    c_r(:,:,ir) = c;
    ell_r(:,:,ir) = ell;
    
    % --- KEY CALCULATION: AGGREGATE ASSET SUPPLY S(r) ---
    S(ir) = g(:,1)'*a*da + g(:,2)'*a*da;
    
    % Compute aggregate effective labor (for reference; add outer loop for full GE)
    L_s(ir) = da * (g(:,1)' * (z(1) * ell(:,1)) + g(:,2)' * (z(2) * ell(:,2)));
    % Note: For full equilibrium, iterate w until L_s ≈ z_ave * mean(ell)
    
    % Display progress
    if mod(ir, 10) == 0
        fprintf('Completed r iteration %d/%d, S = %.4f, L_s = %.4f\n', ir, Ir, S(ir), L_s(ir))
    end
    
end

toc; % End timer

% =========================================================================
% 7. PLOT EQUILIBRIUM
% =========================================================================

Smax = max(S);
amin1 = amin - 0.02;
aaa = linspace(amin1, Smax, Ir);
rrr = linspace(-0.06, 0.06, Ir);
KD_plot = (al*Aprod./(max(rrr + d, 0))).^(1/(1 - al))*z_ave;

figure('Position', [100, 100, 800, 600])
set(gca, 'FontSize', 14)
plot(S, r_grid, 'b-', 'LineWidth', 2.5, 'DisplayName', 'S(r) - Asset Supply')
hold on
plot(KD_plot, rrr, 'r-', 'LineWidth', 2.5, 'DisplayName', 'K^D(r) - Capital Demand')
plot(zeros(1,Ir) + amin, rrr, 'k--', 'LineWidth', 1.5)
plot(aaa, ones(1,Ir)*rho, 'k--', 'LineWidth', 1.5)
plot(aaa, ones(1,Ir)*(-d), 'k--', 'LineWidth', 1.5)

text(0.05, 0.052, '$r = \rho$', 'FontSize', 16, 'interpreter', 'latex')
text(0.05, -0.054, '$r = -\delta$', 'FontSize', 16, 'interpreter', 'latex')
text(0.1, 0.02, '$S(r)$', 'FontSize', 16, 'interpreter', 'latex')
text(0.29, 0.02, '$r=F_K(K,1)-\delta$', 'FontSize', 16, 'interpreter', 'latex')
text(0.01, 0.035, '$a=\underline{a}$', 'FontSize', 16, 'interpreter', 'latex')

ylabel('$r$', 'FontSize', 16, 'interpreter', 'latex')
xlabel('$K$', 'FontSize', 16, 'interpreter', 'latex')
title('Asset Supply and Demand with Endogenous Labor', 'FontSize', 18)
legend('show', 'Location', 'best')
grid on
xlim([amin1 0.6])

% Save figure
print('-dpng', 'aiyagari_endog_labor_equilibrium.png', '-r300')

% =========================================================================
% 8. SAVE RESULTS
% =========================================================================

save('aiyagari_endog_labor_results.mat', ...
     'r_grid', 'S', 'KD', 'w_r', 'V_r', 'g_r', 'c_r', 'ell_r', ...
     'adot', 'a', 'z', 'ga', 'Frisch', 'rho', 'al', 'd', 'L_s')

fprintf('\n=== SIMULATION COMPLETE ===\n')
fprintf('Results saved\n')