%% AIYAGARI MODEL WITH ENDOGENOUS LABOR - 2 FIRMS (FORMAL/INFORMAL)
% Based on Achdou, Han, Lasry, Lions, and Moll (2017)
% Extension with 2 sectors: Formal (uses K+L) and Informal (uses L only)
% Author: Adapted for heterogeneous agents with formal/informal duality
% Date: 2026

clear all; clc; close all;

tic; % Start timer

% =========================================================================
% 1. PARAMETERS
% =========================================================================

% --- Household Parameters ---
ga = 2;           % Relative risk aversion (γ)
rho = 0.05;       % Subjective discount rate (ρ)
Frisch = 0.5;     % Frisch elasticity of labor supply (φ)

% --- Two-Firm Parameters ---
psi_F = 1.0;      % Disutility weight formal sector
psi_I = 0.8;      % Disutility weight informal sector (< psi_F: informal "easier")
theta = 0.5;      % Shock attenuation for informal sector (z_I = theta * z)

% --- Productivity States (z) ---
z1 = 0.2;         % Low productivity
z2 = 2*z1;        % High productivity
z = [z1, z2];     % Productivity vector
la1 = 1;          % Poisson intensity z1 -> z2
la2 = 1;          % Poisson intensity z2 -> z1
la = [la1, la2];

% Average productivity (for normalization)
z_ave = (z1*la2 + z2*la1)/(la1 + la2);

% --- Formal Firm Parameters ---
Aprod = 0.3;      % Total factor productivity formal (A_F)
al = 1/3;         % Capital share (α)
d = 0.05;         % Depreciation rate (δ)

% --- Informal Firm Parameters ---
A_I = 0.15;       % Informal sector productivity (w_I = A_I, constant)

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
ell_F = zeros(I, 2); % Formal labor supply
ell_I = zeros(I, 2); % Informal labor supply

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

% --- Formal Firm Problem (given r) ---
KD = (al*Aprod/(r + d))^(1/(1 - al))*z_ave;
w_F = (1 - al)*Aprod*(KD/z_ave)^al;

% --- Informal Firm ---
w_I = A_I;  % Constant marginal product of labor

% --- Initial Value Function Guess ---
v0(:,1) = (w_F*z(1) + w_I*theta*z(1) + max(r,0.01)*a).^(1 - ga)/(1 - ga)/rho;
v0(:,2) = (w_F*z(2) + w_I*theta*z(2) + max(r,0.01)*a).^(1 - ga)/(1 - ga)/rho;

% =========================================================================
% 4. OUTER LOOP: ITERATE OVER INTEREST RATES
% =========================================================================

fprintf('Starting 2-Firm Aiyagari Model...\n');
fprintf('Parameters: theta=%.2f, A_I=%.2f, psi_I=%.2f\n', theta, A_I, psi_I);

for ir = 1:Ir

    r = r_grid(ir);

    % --- 4.1. Firm Prices (given r) ---
    % Formal firm
    KD(ir) = (al*Aprod/(r + d))^(1/(1 - al))*z_ave;
    w_F = (1 - al)*Aprod*(KD(ir)/z_ave)^al;
    w_F_r(ir) = w_F;

    % Informal firm (constant wage)
    w_I = A_I;
    w_I_r(ir) = w_I;

    % --- 4.2. Warm Start ---
    if ir > 1
        v0 = V_r(:,:,ir - 1);
    end

    v = v0;

    % =====================================================================
    % 5. INNER LOOP: SOLVE HJB (given r, w_F, w_I)
    % =====================================================================

    for n = 1:maxit
        V = v;

        % --- 5.1. Finite Difference Approximations ---
        % Forward derivative
        dVf(1:I-1,:) = (V(2:I,:) - V(1:I-1,:))/da;
        % Boundary at amax: state constraint (consumption equals income)
        dVf(I,:) = (w_F*z + w_I*theta*z + r*amax).^(-ga);

        % Backward derivative
        dVb(2:I,:) = (V(2:I,:) - V(1:I-1,:))/da;
        % Boundary at amin: state constraint (no borrowing)
        dVb(1,:) = (w_F*z + w_I*theta*z + r*amin).^(-ga);

        % Ensure positive derivatives (CRITICAL: prevents complex numbers in FOCs)
        dVf = max(dVf, 1e-10);
        dVb = max(dVb, 1e-10);

        % --- 5.2. Optimal Policies from FOCs ---

        % Forward difference policies
        cf = dVf.^(-1/ga);
        ell_Ff = ((dVf .* w_F .* zz) / psi_F).^Frisch;
        ell_If = ((dVf .* w_I .* theta .* zz) / psi_I).^Frisch;
        ssf = w_F.*zz.*ell_Ff + w_I.*theta.*zz.*ell_If + r.*aa - cf;

        % Backward difference policies
        cb = dVb.^(-1/ga);
        ell_Fb = ((dVb .* w_F .* zz) / psi_F).^Frisch;
        ell_Ib = ((dVb .* w_I .* theta .* zz) / psi_I).^Frisch;
        ssb = w_F.*zz.*ell_Fb + w_I.*theta.*zz.*ell_Ib + r.*aa - cb;

        % Zero drift policies (consume all income)
        c0 = w_F.*zz.*ell_Ff + w_I.*theta.*zz.*ell_If + r.*aa;
        c0 = max(c0, 1e-10);
        ell_F0 = ell_Ff;  % Use forward as approximation
        ell_I0 = ell_If;

        % --- 5.3. Upwind Scheme ---
        If = ssf > 0;
        Ib = ssb < 0 & ~If;
        I0 = ~(If | Ib);

        % Upwind policies
        c = cf.*If + cb.*Ib + c0.*I0;
        ell_F = ell_Ff.*If + ell_Fb.*Ib + ell_F0.*I0;
        ell_I = ell_If.*If + ell_Ib.*Ib + ell_I0.*I0;

        % Ensure positive consumption
        c = max(c, 1e-10);

        % Utility with separate disutility
        u = c.^(1-ga)/(1-ga) ...
            - psi_F * ell_F.^(1+1/Frisch)/(1+1/Frisch) ...
            - psi_I * ell_I.^(1+1/Frisch)/(1+1/Frisch);

        % --- 5.4. Transition Matrix A (using min/max formulation) ---
        % This formulation is mathematically equivalent to Ib*ssb but more robust
        % See aiyagari_poisson_steadystate.m lines 156-158

        % Compute upwind savings for matrix construction
        ss_Upwind = ssf.*If + ssb.*Ib;  % Combined drift

        % Coefficients for tridiagonal matrix
        % X: coefficient for v_{i-1} (comes from backward when saving < 0)
        X = -min(ss_Upwind, 0) / da;
        % Y: diagonal coefficient
        Y = -max(ss_Upwind, 0) / da + min(ss_Upwind, 0) / da;
        % Z: coefficient for v_{i+1} (comes from forward when saving > 0)
        Z = max(ss_Upwind, 0) / da;

        % Build sparse matrices for each productivity state
        % Note: spdiags with [X(2:I,1); 0] pads with zero at end for -1 diagonal
        A1 = spdiags(Y(:,1), 0, I, I) + ...
            spdiags([X(2:I,1); 0], -1, I, I) + ...
            spdiags([0; Z(1:I-1,1)], 1, I, I);
        A2 = spdiags(Y(:,2), 0, I, I) + ...
            spdiags([X(2:I,2); 0], -1, I, I) + ...
            spdiags([0; Z(1:I-1,2)], 1, I, I);

        % Full transition matrix including Poisson jumps between states
        A = [A1, sparse(I,I); sparse(I,I), A2] + Aswitch;

        % Verify and adjust row sums (generator rows must sum to 0)
        row_sum = sum(A, 2);
        if max(abs(row_sum)) > 1e-9
            % Adjust diagonal to enforce row sums = 0
            for j = 1:2*I
                A(j,j) = A(j,j) - row_sum(j);
            end
            % Verify adjustment worked
            if max(abs(sum(A,2))) > 1e-9
                fprintf('Warning: Row sum adjustment failed at r=%.4f, max error=%.2e\n', r, max(abs(sum(A,2))));
            end
        end

        % --- 5.5. Solve Implicit System ---
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

    % Ensure non-negative distribution
    gg = max(gg, 0);

    % Normalize (total mass = 1)
    g_sum = gg'*ones(2*I,1)*da;
    if g_sum > 1e-10
        gg = gg/g_sum;
    else
        % Fallback: uniform distribution
        gg = ones(2*I,1)/(2*I*da);
    end

    g = [gg(1:I), gg(I+1:2*I)];

    % --- 6.1. Store Results ---
    g_r(:,:,ir) = g;
    V_r(:,:,ir) = V;
    c_r(:,:,ir) = c;
    ell_F_r(:,:,ir) = ell_F;
    ell_I_r(:,:,ir) = ell_I;

    % --- 6.2. Aggregate Variables ---
    % Asset supply
    S(ir) = g(:,1)'*a*da + g(:,2)'*a*da;

    % Aggregate labor formal
    L_F(ir) = da * (g(:,1)' * (z(1) * ell_F(:,1)) + g(:,2)' * (z(2) * ell_F(:,2)));

    % Aggregate labor informal (with theta attenuation)
    L_I(ir) = da * (g(:,1)' * (theta * z(1) * ell_I(:,1)) + g(:,2)' * (theta * z(2) * ell_I(:,2)));

    % Display progress
    if mod(ir, 10) == 0
        fprintf('Completed r iteration %d/%d, S = %.4f, L_F = %.4f, L_I = %.4f\n', ...
            ir, Ir, S(ir), L_F(ir), L_I(ir))
    end

end

toc; % End timer

% =========================================================================
% 7. PLOT RESULTS
% =========================================================================

% Create output directory
output_dir = 'output_graphs';
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

% --- Figure 1: Asset Supply and Demand ---
figure('Position', [100, 100, 800, 600])
Smax = max(S);
amin1 = amin - 0.02;
aaa = linspace(amin1, Smax, Ir);
rrr = linspace(-0.06, 0.06, Ir);
KD_plot = (al*Aprod./(max(rrr + d, 0))).^(1/(1 - al))*z_ave;

plot(S, r_grid, 'b-', 'LineWidth', 2.5)
hold on
plot(KD_plot, rrr, 'r-', 'LineWidth', 2.5)
plot(zeros(1,Ir) + amin, rrr, 'k--', 'LineWidth', 1.5)
plot(aaa, ones(1,Ir)*rho, 'k--', 'LineWidth', 1.5)
plot(aaa, ones(1,Ir)*(-d), 'k--', 'LineWidth', 1.5)

ylabel('$r$', 'FontSize', 16, 'interpreter', 'latex')
xlabel('$K$', 'FontSize', 16, 'interpreter', 'latex')
title('Asset Supply and Demand - 2 Firms', 'FontSize', 18)
legend('S(r) - Asset Supply', 'K^D(r) - Capital Demand', 'Location', 'best')
grid on
xlim([amin1 0.6])
set(gca, 'FontSize', 14)
print('-dpng', [output_dir '/equilibrium_2firms.png'], '-r300')

% --- Figure 2: Labor Supply by Sector ---
figure('Position', [100, 100, 800, 600])
subplot(2,1,1)
plot(r_grid, L_F, 'b-', 'LineWidth', 2)
hold on
plot(r_grid, L_I, 'r--', 'LineWidth', 2)
xlabel('$r$', 'FontSize', 14, 'interpreter', 'latex')
ylabel('Labor Supply', 'FontSize', 14)
title('Aggregate Labor by Sector', 'FontSize', 16)
legend('L_F (Formal)', 'L_I (Informal)', 'Location', 'best')
grid on
set(gca, 'FontSize', 12)

subplot(2,1,2)
plot(r_grid, L_I./(L_F + L_I), 'k-', 'LineWidth', 2)
xlabel('$r$', 'FontSize', 14, 'interpreter', 'latex')
ylabel('Informal Share', 'FontSize', 14)
title('Informal Labor Share: $L_I / (L_F + L_I)$', 'FontSize', 16, 'interpreter', 'latex')
grid on
set(gca, 'FontSize', 12)
print('-dpng', [output_dir '/labor_by_sector.png'], '-r300')

% --- Figure 3: Labor Policy Functions (for middle r) ---
ir_mid = round(Ir/2);
figure('Position', [100, 100, 1000, 400])
subplot(1,2,1)
plot(a, ell_F_r(:,1,ir_mid), 'b-', 'LineWidth', 2)
hold on
plot(a, ell_F_r(:,2,ir_mid), 'b--', 'LineWidth', 2)
xlabel('Wealth $a$', 'FontSize', 14, 'interpreter', 'latex')
ylabel('$\ell_F(a,z)$', 'FontSize', 14, 'interpreter', 'latex')
title('Formal Labor Supply', 'FontSize', 16)
legend('z = z_1 (low)', 'z = z_2 (high)', 'Location', 'best')
grid on
set(gca, 'FontSize', 12)

subplot(1,2,2)
plot(a, ell_I_r(:,1,ir_mid), 'r-', 'LineWidth', 2)
hold on
plot(a, ell_I_r(:,2,ir_mid), 'r--', 'LineWidth', 2)
xlabel('Wealth $a$', 'FontSize', 14, 'interpreter', 'latex')
ylabel('$\ell_I(a,z)$', 'FontSize', 14, 'interpreter', 'latex')
title('Informal Labor Supply', 'FontSize', 16)
legend('z = z_1 (low)', 'z = z_2 (high)', 'Location', 'best')
grid on
set(gca, 'FontSize', 12)
print('-dpng', [output_dir '/labor_policies.png'], '-r300')

% --- Figure 4: Informal Share by Wealth ---
figure('Position', [100, 100, 800, 500])
informal_share = ell_I_r(:,:,ir_mid) ./ (ell_F_r(:,:,ir_mid) + ell_I_r(:,:,ir_mid));
plot(a, informal_share(:,1), 'k-', 'LineWidth', 2)
hold on
plot(a, informal_share(:,2), 'k--', 'LineWidth', 2)
xlabel('Wealth $a$', 'FontSize', 14, 'interpreter', 'latex')
ylabel('Informal Share $\ell_I / (\ell_F + \ell_I)$', 'FontSize', 14, 'interpreter', 'latex')
title('Informal Labor Share by Wealth', 'FontSize', 16)
legend('z = z_1 (low)', 'z = z_2 (high)', 'Location', 'best')
grid on
set(gca, 'FontSize', 12)
print('-dpng', [output_dir '/informal_share_by_wealth.png'], '-r300')

% =========================================================================
% 8. SAVE RESULTS
% =========================================================================

save('aiyagari_2firms_results.mat', ...
    'r_grid', 'S', 'KD', 'w_F_r', 'w_I_r', 'V_r', 'g_r', 'c_r', ...
    'ell_F_r', 'ell_I_r', 'L_F', 'L_I', ...
    'a', 'z', 'ga', 'Frisch', 'rho', 'al', 'd', ...
    'theta', 'A_I', 'psi_F', 'psi_I')

fprintf('\n=== SIMULATION COMPLETE ===\n')
fprintf('Results saved to aiyagari_2firms_results.mat\n')
fprintf('Graphs saved to %s/\n', output_dir)
