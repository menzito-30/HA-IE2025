%% AIYAGARI MODEL - 2 FIRMS - VERSION 5: MIT SHOCK
% Transition dynamics after unexpected productivity shock
% Based on Moll's aiyagari_poisson_MITshock.m adapted for 2 firms
%
% This script:
%   1. First runs v4 in EQUILIBRIUM_MODE=2 to get steady state
%   2. Then simulates MIT shock transition dynamics
%
% Author: Macroeconomic HA expert
% Date: 2026

clear all; clc; close all;

tic;

% =========================================================================
% 1. PARAMETERS (Same as v4)
% =========================================================================

% Household
ga = 2;
rho = 0.05;
Frisch = 0.5;

% Two-Firm
psi_F = 1.0;
psi_I = 0.8;
theta = 0.5;

% Productivity states
z1 = 0.2;
z2 = 2*z1;
z = [z1, z2];
la1 = 1; la2 = 1;
la = [la1, la2];
z_ave = (z1*la2 + z2*la1)/(la1 + la2);

% Formal Firm
Aprod = 0.3;
al = 1/3;
d = 0.05;

% Informal Firm
A_I = 0.15;

% =========================================================================
% 2. MIT SHOCK PARAMETERS (CONFIGURABLE)
% =========================================================================

% Time parameters
T = 200;              % Time horizon (years)
N = 400;              % Number of time steps
dt = T/N;             % Time step
time = (0:N-1)*dt;    % Time grid

% -----------------------------------------------------------------
% SHOCK CONFIGURATION (MODIFY THESE)
% -----------------------------------------------------------------

% Shock sign: negative = recession, positive = boom
shock_formal = -0.03;    % -3% drop in formal productivity (like Moll)
shock_informal = -0.03;  % -3% drop in informal productivity
% Set to 0 for no shock to that sector

% Persistence
corr = 0.8;              % Persistence (higher = slower recovery)
nu = 1 - corr;           % Mean reversion speed

% -----------------------------------------------------------------

% Construct TFP sequences
A_F_st = Aprod;          % Steady state formal productivity
A_I_st = A_I;            % Steady state informal productivity

% Formal sector shock
A_F_t = zeros(N,1);
A_F_t(1) = (1 + shock_formal) * A_F_st;
for n = 1:N-1
    A_F_t(n+1) = dt * nu * (A_F_st - A_F_t(n)) + A_F_t(n);
end

% Informal sector shock
A_I_t = zeros(N,1);
A_I_t(1) = (1 + shock_informal) * A_I_st;
for n = 1:N-1
    A_I_t(n+1) = dt * nu * (A_I_st - A_I_t(n)) + A_I_t(n);
end

fprintf('=== MIT SHOCK CONFIGURATION ===\n');
fprintf('Formal shock:   %.1f%%\n', shock_formal*100);
fprintf('Informal shock: %.1f%%\n', shock_informal*100);
fprintf('Persistence:    %.2f\n', corr);
fprintf('===============================\n\n');

% Price iteration parameters
max_price_it = 2000;          % Increased for better convergence
convergence_criterion = 1e-3; % Slightly relaxed
relax_K = 0.05;               % Relaxation for capital
relax_L = 0.005;              % VERY conservative for labor (avoid instability)

% =========================================================================
% 3. GRIDS
% =========================================================================

I = 500;              % Increased for smoother inequality measures
amin = 0;
amax = 20;
a = linspace(amin, amax, I)';
da = (amax - amin)/(I - 1);

aa = [a, a];
zz = ones(I,1)*z;

maxit = 100;
crit = 10^(-6);
Delta = 1000;

Aswitch = [-speye(I)*la(1), speye(I)*la(1);
    speye(I)*la(2), -speye(I)*la(2)];

% =========================================================================
% PHASE 1: STEADY STATE (Using v4 approach)
% =========================================================================

fprintf('=== PHASE 1: SOLVING STEADY STATE ===\n\n');

% Bisection for steady state
r_low = -0.04;
r_high = 0.045;
tol_r = 1e-5;
max_bisect = 50;

r = (r_low + r_high)/2;
KD_init = (al*Aprod/(r + d))^(1/(1 - al))*z_ave;
w_F = (1 - al)*Aprod*(KD_init/z_ave)^al;
w_I = A_I;

v0(:,1) = (w_F*z(1) + w_I*theta*z(1) + max(r,0.01)*a).^(1-ga)/(1-ga)/rho;
v0(:,2) = (w_F*z(2) + w_I*theta*z(2) + max(r,0.01)*a).^(1-ga)/(1-ga)/rho;

for iter = 1:max_bisect
    r = (r_low + r_high)/2;

    % Call robust solver
    [S_mid, KD_mid, w_F_mid, L_F_mid, L_I_mid, V, g, c, ell_F, ell_I, A_ss, v0] = ...
        solve_steady_state_robust(r, v0, a, z, la, ga, rho, Frisch, psi_F, psi_I, theta, ...
        Aprod, al, d, A_I, z_ave, I, da, aa, zz, maxit, crit, Delta, Aswitch);

    excess = S_mid - KD_mid;

    if mod(iter, 5) == 0
        fprintf('Bisect iter %2d: r=%.6f, S=%.4f, K^D=%.4f, excess=%.4f\n', ...
            iter, r, S_mid, KD_mid, excess);
    end

    if abs(excess) < tol_r || (r_high - r_low) < tol_r
        fprintf('*** STEADY STATE FOUND ***\n');
        break;
    end

    if excess > 0
        r_high = r;
    else
        r_low = r;
    end
end

% Save steady state
r_st = r;
K_st = KD_mid;
w_F_st = w_F_mid;
w_I_st = A_I;
L_F_st = L_F_mid;
L_I_st = L_I_mid;
V_st = V;
gg_st = [g(:,1); g(:,2)];
c_st = c;
ell_F_st = ell_F;
ell_I_st = ell_I;
A_st = A_ss;

fprintf('r* = %.4f, K* = %.4f, L_F* = %.4f, L_I* = %.4f\n\n', r_st, K_st, L_F_st, L_I_st);
toc;

% =========================================================================
% PHASE 2: TRANSITION DYNAMICS
% =========================================================================

fprintf('\n=== PHASE 2: TRANSITION DYNAMICS ===\n\n');

% Preallocation
gg = cell(N+1, 1);
K_t = zeros(N, 1);
K_out = zeros(N, 1);
L_F_out = zeros(N, 1);  % Output labor aggregates
L_I_out = zeros(N, 1);
r_t = zeros(N, 1);
w_F_t_vec = zeros(N, 1);
w_I_t_vec = zeros(N, 1);
L_F_t_vec = zeros(N, 1);
L_I_t_vec = zeros(N, 1);
C_t = zeros(N, 1);
dist_it = zeros(max_price_it, 1);

% Initial guess: capital and labor - IMPROVED INITIALIZATION
% Use linear interpolation from shocked state to steady state
% This provides a better starting point than constant K_st
K_shock_impact = K_st * (1 + shock_formal);  % Approximate impact of shock on K
for n = 1:N
    % Smooth transition from K_shock_impact to K_st
    weight = 1 - exp(-nu * time(n));  % Same dynamics as TFP recovery
    K_t(n) = K_shock_impact + (K_st - K_shock_impact) * weight;
end
L_F_t_vec = L_F_st * ones(N, 1);  % Initialize with steady state labor
L_I_t_vec = L_I_st * ones(N, 1);

dVf = zeros(I, 2);
dVb = zeros(I, 2);

% Price iteration (over K and L together)
for it = 1:max_price_it
    fprintf('Price Iteration %d', it);

    % Compute prices from capital AND labor path (with shock)
    % CORRECTED: Use L_F_t instead of z_ave
    for n = 1:N
        % Formal firm: F(K,L) = A_F * K^alpha * L_F^(1-alpha)
        L_F_n = max(L_F_t_vec(n), 0.01);  % Avoid division by zero
        w_F_t_vec(n) = (1 - al) * A_F_t(n) * (K_t(n)/L_F_n)^al;
        r_t(n) = al * A_F_t(n) * (K_t(n)/L_F_n)^(al-1) - d;

        % Informal firm: w_I = A_I (linear technology)
        w_I_t_vec(n) = A_I_t(n);
    end

    % -----------------------------------------------------------------
    % BACKWARD HJB: from t=T to t=0
    % -----------------------------------------------------------------
    V = V_st;
    A_t = cell(N, 1);
    c_t = zeros(I, 2, N);
    ell_F_t_mat = zeros(I, 2, N);
    ell_I_t_mat = zeros(I, 2, N);

    for n = N:-1:1
        w_F = w_F_t_vec(n);
        w_I = w_I_t_vec(n);
        r = r_t(n);

        % Finite differences
        dVf(1:I-1,:) = (V(2:I,:) - V(1:I-1,:))/da;
        dVb(2:I,:) = (V(2:I,:) - V(1:I-1,:))/da;

        % Boundaries
        for j = 1:2
            income_max = w_F*z(j) + w_I*theta*z(j) + r*amax;
            dVf(I,j) = max(income_max, 1e-10)^(-ga);
            income_min = w_F*z(j) + w_I*theta*z(j) + r*amin;
            dVb(1,j) = max(income_min, 1e-10)^(-ga);
        end

        % Policies from FOCs
        dVf_pos = max(dVf, 1e-10);
        dVb_pos = max(dVb, 1e-10);

        cf = dVf_pos.^(-1/ga);
        cb = dVb_pos.^(-1/ga);

        % Labor from FOCs
        ell_Ff = ((dVf_pos .* w_F .* zz) / psi_F).^Frisch;
        ell_If = ((dVf_pos .* w_I * theta .* zz) / psi_I).^Frisch;
        ell_Fb = ((dVb_pos .* w_F .* zz) / psi_F).^Frisch;
        ell_Ib = ((dVb_pos .* w_I * theta .* zz) / psi_I).^Frisch;

        ssf = w_F*zz.*ell_Ff + w_I*theta*zz.*ell_If + r*aa - cf;
        ssb = w_F*zz.*ell_Fb + w_I*theta*zz.*ell_Ib + r*aa - cb;

        % Zero drift
        c0 = w_F*zz + w_I*theta*zz + r*aa;
        c0 = max(c0, 1e-10);
        dV0 = c0.^(-ga);
        ell_F0 = ((dV0 .* w_F .* zz) / psi_F).^Frisch;
        ell_I0 = ((dV0 .* w_I * theta .* zz) / psi_I).^Frisch;

        % Upwind
        If = ssf > 0;
        Ib = ssb < 0 & ~If;
        I0 = ~(If | Ib);

        c = cf.*If + cb.*Ib + c0.*I0;
        ell_F_n = ell_Ff.*If + ell_Fb.*Ib + ell_F0.*I0;
        ell_I_n = ell_If.*If + ell_Ib.*Ib + ell_I0.*I0;

        c_t(:,:,n) = c;
        ell_F_t_mat(:,:,n) = ell_F_n;
        ell_I_t_mat(:,:,n) = ell_I_n;

        u = c.^(1-ga)/(1-ga) - psi_F*ell_F_n.^(1+1/Frisch)/(1+1/Frisch) ...
            - psi_I*ell_I_n.^(1+1/Frisch)/(1+1/Frisch);

        % Transition matrix
        ss_Upwind = ssf.*If + ssb.*Ib;
        X = -min(ss_Upwind, 0) / da;
        Y = -max(ss_Upwind, 0) / da + min(ss_Upwind, 0) / da;
        Z = max(ss_Upwind, 0) / da;

        A1 = spdiags(Y(:,1), 0, I, I) + spdiags([X(2:I,1); 0], -1, I, I) + ...
            spdiags([0; Z(1:I-1,1)], 1, I, I);
        A2 = spdiags(Y(:,2), 0, I, I) + spdiags([X(2:I,2); 0], -1, I, I) + ...
            spdiags([0; Z(1:I-1,2)], 1, I, I);
        A = [A1, sparse(I,I); sparse(I,I), A2] + Aswitch;

        A_t{n} = A;

        % Implicit update
        B = (1/dt + rho)*speye(2*I) - A;
        u_stacked = [u(:,1); u(:,2)];
        V_stacked = [V(:,1); V(:,2)];
        b = u_stacked + V_stacked/dt;
        V_stacked = B\b;
        V = [V_stacked(1:I), V_stacked(I+1:2*I)];
    end

    % -----------------------------------------------------------------
    % FORWARD KFE: from t=0 to t=T
    % -----------------------------------------------------------------
    gg{1} = gg_st;

    for n = 1:N
        AT = A_t{n}';
        gg{n+1} = (speye(2*I) - AT*dt) \ gg{n};

        % Normalize
        g_sum = sum(gg{n+1}) * da;
        if g_sum > 1e-10
            gg{n+1} = gg{n+1} / g_sum;
        end

        % Aggregate capital
        K_out(n) = gg{n}(1:I)'*a*da + gg{n}(I+1:2*I)'*a*da;

        % Aggregate consumption
        c_n = c_t(:,:,n);
        C_t(n) = gg{n}' * [c_n(:,1); c_n(:,2)] * da;

        % Aggregate labor (compute output values)
        ell_F_n = ell_F_t_mat(:,:,n);
        ell_I_n = ell_I_t_mat(:,:,n);
        L_F_out(n) = da * (gg{n}(1:I)' * (z(1)*ell_F_n(:,1)) + gg{n}(I+1:2*I)' * (z(2)*ell_F_n(:,2)));
        L_I_out(n) = da * (gg{n}(1:I)' * (theta*z(1)*ell_I_n(:,1)) + gg{n}(I+1:2*I)' * (theta*z(2)*ell_I_n(:,2)));
    end

    % Check convergence (on K and L)
    dist_K = max(abs(K_out - K_t));
    dist_L_F = max(abs(L_F_out - L_F_t_vec));
    dist_L_I = max(abs(L_I_out - L_I_t_vec));
    dist_it(it) = dist_K;  % Track K convergence
    fprintf('  K_dist=%.4f, L_F_dist=%.4f, L_I_dist=%.4f\n', dist_K, dist_L_F, dist_L_I);

    if dist_K < convergence_criterion && dist_L_F < convergence_criterion
        fprintf('\n*** TRANSITION DYNAMICS CONVERGED ***\n');
        break;
    end

    % Update capital AND labor paths with SEPARATE relaxation
    K_t = relax_K * K_out + (1 - relax_K) * K_t;
    L_F_t_vec = relax_L * L_F_out + (1 - relax_L) * L_F_t_vec;
    L_I_t_vec = relax_L * L_I_out + (1 - relax_L) * L_I_t_vec;
end

toc;

% =========================================================================
% PHASE 3: ANALYSIS AND PLOTS
% =========================================================================

fprintf('\n=== PHASE 3: GENERATING PLOTS ===\n\n');

output_dir = 'output_graphs';
if ~exist(output_dir, 'dir'), mkdir(output_dir); end

% Calculate informal share over time
informal_share_t = L_I_t_vec ./ (L_F_t_vec + L_I_t_vec + 1e-10);

% =========================================================================
% INEQUALITY MEASURES (following Moll)
% =========================================================================

fprintf('Calculating inequality measures...\n');

Wealth_Gini_t = zeros(N,1);
Income_Gini_t = zeros(N,1);
Formal_Inc_Gini_t = zeros(N,1);
top_wealth_10_t = zeros(N,1);
top_income_10_t = zeros(N,1);

for n = 1:N
    % Get distribution at time n
    g_n = gg{n};

    % Marginal wealth distribution (sum over z states)
    g_a_cont = g_n(1:I) + g_n(I+1:2*I);
    g_a = g_a_cont * da;

    % Wealth Gini
    S_a = cumsum(g_a .* a) / max(sum(g_a .* a), 1e-10);
    trapez_a = (1/2) * (S_a(1)*g_a(1) + sum((S_a(2:I) + S_a(1:I-1)) .* g_a(2:I)));
    Wealth_Gini_t(n) = 1 - 2*trapez_a;

    % Top 10% Wealth Share
    G_a = cumsum(g_a);
    [~, idx] = min(abs((1-G_a) - 0.1));
    top_wealth_10_t(n) = 1 - S_a(idx);

    % Total income: formal + informal + asset income
    % y(a,z) = w_F*z*ell_F(a,z) + w_I*theta*z*ell_I(a,z) + r*a
    ell_F_n = ell_F_t_mat(:,:,n);
    ell_I_n = ell_I_t_mat(:,:,n);
    w_F = w_F_t_vec(n);
    w_I = w_I_t_vec(n);
    r = r_t(n);

    % Income for each (a,z) pair
    y_az = zeros(I, 2);
    for j = 1:2
        y_az(:,j) = w_F * z(j) * ell_F_n(:,j) + w_I * theta * z(j) * ell_I_n(:,j) + r * a;
    end

    % Sort by income
    Ny = 2*I;
    yy = reshape(y_az, Ny, 1);
    [yy_sorted, index] = sort(yy);
    g_y = g_n(index) * da;

    % Income Gini
    S_y = cumsum(g_y .* yy_sorted) / max(sum(g_y .* yy_sorted), 1e-10);
    trapez_y = (1/2) * (S_y(1)*g_y(1) + sum((S_y(2:Ny) + S_y(1:Ny-1)) .* g_y(2:Ny)));
    Income_Gini_t(n) = 1 - 2*trapez_y;

    % Top 10% Income Share
    G_y = cumsum(g_y);
    [~, idx] = min(abs((1-G_y) - 0.1));
    top_income_10_t(n) = 1 - S_y(idx);

    % Formal income only Gini
    y_formal = zeros(I, 2);
    for j = 1:2
        y_formal(:,j) = w_F * z(j) * ell_F_n(:,j);
    end
    yf = reshape(y_formal, Ny, 1);
    [yf_sorted, index_f] = sort(yf);
    g_yf = g_n(index_f) * da;
    S_yf = cumsum(g_yf .* yf_sorted) / max(sum(g_yf .* yf_sorted), 1e-10);
    trapez_yf = (1/2) * (S_yf(1)*g_yf(1) + sum((S_yf(2:Ny) + S_yf(1:Ny-1)) .* g_yf(2:Ny)));
    Formal_Inc_Gini_t(n) = 1 - 2*trapez_yf;
end

fprintf('Inequality measures calculated.\n');

% Figure 1: Shock and main variables
figure('Position', [100, 100, 1200, 800])

subplot(2,3,1)
plot(time, A_F_t, 'b-', time, A_I_t, 'r--', 'LineWidth', 2)
hold on
plot(time, A_F_st*ones(N,1), 'b:', 'LineWidth', 1)
plot(time, A_I_st*ones(N,1), 'r:', 'LineWidth', 1)
xlabel('Year'); ylabel('Productivity')
title('MIT Shock: Productivity')
legend('A_F', 'A_I', 'Location', 'best')
xlim([0 50])
grid on

subplot(2,3,2)
plot(time, K_t, 'b-', 'LineWidth', 2)
hold on
plot(time, K_st*ones(N,1), 'k--', 'LineWidth', 1)
xlabel('Year'); ylabel('K')
title('Aggregate Capital')
xlim([0 50])
grid on

subplot(2,3,3)
plot(time, w_F_t_vec, 'b-', time, w_I_t_vec, 'r--', 'LineWidth', 2)
hold on
plot(time, w_F_st*ones(N,1), 'b:', 'LineWidth', 1)
plot(time, w_I_st*ones(N,1), 'r:', 'LineWidth', 1)
xlabel('Year'); ylabel('Wage')
title('Wages')
legend('w_F', 'w_I', 'Location', 'best')
xlim([0 50])
grid on

subplot(2,3,4)
plot(time, r_t, 'b-', 'LineWidth', 2)
hold on
plot(time, r_st*ones(N,1), 'k--', 'LineWidth', 1)
xlabel('Year'); ylabel('r')
title('Interest Rate')
xlim([0 50])
grid on

subplot(2,3,5)
plot(time, L_F_t_vec, 'b-', time, L_I_t_vec, 'r--', 'LineWidth', 2)
hold on
plot(time, L_F_st*ones(N,1), 'b:', 'LineWidth', 1)
plot(time, L_I_st*ones(N,1), 'r:', 'LineWidth', 1)
xlabel('Year'); ylabel('Labor')
title('Formal vs Informal Labor')
legend('L_F', 'L_I', 'Location', 'best')
xlim([0 50])
grid on

subplot(2,3,6)
plot(time, 100*informal_share_t, 'k-', 'LineWidth', 2)
hold on
plot(time, 100*L_I_st/(L_F_st+L_I_st)*ones(N,1), 'k--', 'LineWidth', 1)
xlabel('Year'); ylabel('Share (%)')
title('Informal Share of Labor')
xlim([0 50])
grid on

print('-dpng', [output_dir '/mit_shock_dynamics.png'], '-r300')

% Figure 2: Convergence
figure('Position', [100, 100, 600, 400])
plot(1:it, dist_it(1:it), 'b-o', 'LineWidth', 2)
xlabel('Iteration'); ylabel('Max |K_{out} - K_t|')
title('Price Iteration Convergence')
grid on
print('-dpng', [output_dir '/mit_shock_convergence.png'], '-r300')

% Figure 3: Inequality Measures (like Moll)
figure('Position', [100, 100, 1200, 600])

subplot(2,3,1)
plot(time, Wealth_Gini_t, 'b-', 'LineWidth', 2)
hold on
plot(time, Wealth_Gini_t(end)*ones(N,1), 'k--', 'LineWidth', 1)
xlabel('Year'); ylabel('Gini')
title('Wealth Gini')
xlim([0 100])
grid on

subplot(2,3,2)
plot(time, Income_Gini_t, 'b-', 'LineWidth', 2)
hold on
plot(time, Income_Gini_t(end)*ones(N,1), 'k--', 'LineWidth', 1)
xlabel('Year'); ylabel('Gini')
title('Total Income Gini')
xlim([0 100])
grid on

subplot(2,3,3)
plot(time, Formal_Inc_Gini_t, 'b-', 'LineWidth', 2)
hold on
plot(time, Formal_Inc_Gini_t(end)*ones(N,1), 'k--', 'LineWidth', 1)
xlabel('Year'); ylabel('Gini')
title('Formal Income Gini')
xlim([0 100])
grid on

subplot(2,3,4)
plot(time, 100*top_wealth_10_t, 'b-', 'LineWidth', 2)
hold on
plot(time, 100*top_wealth_10_t(end)*ones(N,1), 'k--', 'LineWidth', 1)
xlabel('Year'); ylabel('Share (%)')
title('Top 10% Wealth Share')
xlim([0 100])
grid on

subplot(2,3,5)
plot(time, 100*top_income_10_t, 'b-', 'LineWidth', 2)
hold on
plot(time, 100*top_income_10_t(end)*ones(N,1), 'k--', 'LineWidth', 1)
xlabel('Year'); ylabel('Share (%)')
title('Top 10% Income Share')
xlim([0 100])
grid on

subplot(2,3,6)
plot(time, 100*informal_share_t, 'b-', 'LineWidth', 2)
hold on
plot(time, 100*L_I_st/(L_F_st+L_I_st)*ones(N,1), 'k--', 'LineWidth', 1)
xlabel('Year'); ylabel('Share (%)')
title('Informal Labor Share')
xlim([0 100])
grid on

print('-dpng', [output_dir '/mit_shock_inequality.png'], '-r300')

% =========================================================================
% SAVE RESULTS
% =========================================================================

save('aiyagari_2firms_MITshock.mat', ...
    'time', 'A_F_t', 'A_I_t', 'K_t', 'r_t', 'w_F_t_vec', 'w_I_t_vec', ...
    'L_F_t_vec', 'L_I_t_vec', 'C_t', ...
    'r_st', 'K_st', 'w_F_st', 'w_I_st', 'L_F_st', 'L_I_st', ...
    'Wealth_Gini_t', 'Income_Gini_t', 'Formal_Inc_Gini_t', ...
    'top_wealth_10_t', 'top_income_10_t', 'informal_share_t', ...
    'shock_formal', 'shock_informal', 'corr', 'T', 'N')

fprintf('\n=== MIT SHOCK SIMULATION COMPLETE ===\n')
fprintf('Results saved to aiyagari_2firms_MITshock.mat\n')
fprintf('Graphs saved to %s/\n', output_dir)

toc;


% =========================================================================
% HELPER FUNCTION: Solve steady state for given r (copied from v4)
% =========================================================================

function [S, KD, w_F, L_F, L_I, V, g, c, ell_F, ell_I, A, v0_out] = ...
    solve_steady_state_robust(r, v0, a, z, la, ga, rho, Frisch, psi_F, psi_I, theta, ...
    Aprod, al, d, A_I, z_ave, I, da, aa, zz, maxit, crit, Delta, Aswitch)

% Firm prices
KD = (al*Aprod/(r + d))^(1/(1 - al))*z_ave;
w_F = (1 - al)*Aprod*(KD/z_ave)^al;
w_I = A_I;

v = v0;
dVf = zeros(I, 2);
dVb = zeros(I, 2);
amin = a(1);
amax = a(end);

% Precompute boundary dV_min using lab_solve_2firms
dV_min = zeros(1,2);
for j = 1:2
    params = [amin, z(j), w_F, w_I, theta, r, ga, Frisch, psi_F, psi_I];
    c_guess = w_F*z(j) + w_I*theta*z(j) + max(r*amin, 0.01);
    try
        c_min = fzero(@(c) lab_solve_2firms(c, params), c_guess, optimset('Display','off','TolFun',1e-10));
    catch
        c_min = c_guess;
    end
    dV_min(j) = max(c_min, 1e-10)^(-ga);
end

% HJB iteration
for n = 1:maxit
    V = v;

    dVf(1:I-1,:) = (V(2:I,:) - V(1:I-1,:))/da;
    dVb(2:I,:) = (V(2:I,:) - V(1:I-1,:))/da;

    % Boundaries using lab_solve_2firms
    for j = 1:2
        params = [amax, z(j), w_F, w_I, theta, r, ga, Frisch, psi_F, psi_I];
        c_guess = w_F*z(j) + w_I*theta*z(j) + max(r*amax, 0.01);
        try
            c_up = fzero(@(c) lab_solve_2firms(c, params), c_guess, optimset('Display','off','TolFun',1e-10));
        catch
            c_up = c_guess;
        end
        c_up = max(c_up, 1e-10);
        u_c = c_up^(-ga);
        ell_F_up = ((u_c * w_F * z(j)) / psi_F)^Frisch;
        ell_I_up = ((u_c * w_I * theta * z(j)) / psi_I)^Frisch;
        income_up = w_F*z(j)*ell_F_up + w_I*theta*z(j)*ell_I_up + r*amax;
        dVf(I,j) = max(income_up, 1e-10)^(-ga);
    end
    dVb(1,:) = max((V(2,:) - V(1,:))/da, dV_min);

    % Forward/Backward policies
    dVf_pos = max(dVf, 1e-10);
    dVb_pos = max(dVb, 1e-10);
    cf = real(dVf_pos .^ (-1/ga));
    cb = real(dVb_pos .^ (-1/ga));

    ell_Ff = zeros(I, 2); ell_If = zeros(I, 2);
    ell_Fb = zeros(I, 2); ell_Ib = zeros(I, 2);
    for j = 1:2
        for i = 1:I
            u_c = dVf_pos(i,j);
            ell_Ff(i,j) = ((u_c * w_F * z(j)) / psi_F)^Frisch;
            ell_If(i,j) = ((u_c * w_I * theta * z(j)) / psi_I)^Frisch;
            u_c = dVb_pos(i,j);
            ell_Fb(i,j) = ((u_c * w_F * z(j)) / psi_F)^Frisch;
            ell_Ib(i,j) = ((u_c * w_I * theta * z(j)) / psi_I)^Frisch;
        end
    end

    ssf = w_F*zz.*ell_Ff + w_I*theta*zz.*ell_If + r*aa - cf;
    ssb = w_F*zz.*ell_Fb + w_I*theta*zz.*ell_Ib + r*aa - cb;

    % Zero drift using lab_solve_2firms
    c0 = zeros(I, 2); ell_F0 = zeros(I, 2); ell_I0 = zeros(I, 2);
    for j = 1:2
        for i = 1:I
            params = [a(i), z(j), w_F, w_I, theta, r, ga, Frisch, psi_F, psi_I];
            c_guess = w_F*z(j) + w_I*theta*z(j) + max(r*a(i), 0.01);
            try
                c0(i,j) = fzero(@(c) lab_solve_2firms(c, params), c_guess, optimset('Display','off','TolFun',1e-10));
            catch
                c0(i,j) = c_guess;
            end
            c0(i,j) = max(c0(i,j), 1e-10);
            u_c = c0(i,j)^(-ga);
            ell_F0(i,j) = ((u_c * w_F * z(j)) / psi_F)^Frisch;
            ell_I0(i,j) = ((u_c * w_I * theta * z(j)) / psi_I)^Frisch;
        end
    end
    dV0 = max(c0 .^ (-ga), 1e-10);

    % Upwind
    If = ssf > 0;
    Ib = ssb < 0 & ~If;
    I0 = ~(If | Ib);

    c = cf.*If + cb.*Ib + c0.*I0;
    ell_F = ell_Ff.*If + ell_Fb.*Ib + ell_F0.*I0;
    ell_I = ell_If.*If + ell_Ib.*Ib + ell_I0.*I0;

    u = c.^(1-ga)/(1-ga) - psi_F*ell_F.^(1+1/Frisch)/(1+1/Frisch) - psi_I*ell_I.^(1+1/Frisch)/(1+1/Frisch);

    % Transition matrix
    ss_Upwind = ssf.*If + ssb.*Ib;
    X = -min(ss_Upwind, 0) / da;
    Y = -max(ss_Upwind, 0) / da + min(ss_Upwind, 0) / da;
    Z = max(ss_Upwind, 0) / da;

    A1 = spdiags(Y(:,1), 0, I, I) + spdiags([X(2:I,1); 0], -1, I, I) + spdiags([0; Z(1:I-1,1)], 1, I, I);
    A2 = spdiags(Y(:,2), 0, I, I) + spdiags([X(2:I,2); 0], -1, I, I) + spdiags([0; Z(1:I-1,2)], 1, I, I);
    A = [A1, sparse(I,I); sparse(I,I), A2] + Aswitch;

    row_sum = sum(A, 2);
    if max(abs(row_sum)) > 1e-9
        Y(:,1) = Y(:,1) - row_sum(1:I) / I;
        Y(:,2) = Y(:,2) - row_sum(I+1:2*I) / I;
        A1 = spdiags(Y(:,1), 0, I, I) + spdiags([X(2:I,1); 0], -1, I, I) + spdiags([0; Z(1:I-1,1)], 1, I, I);
        A2 = spdiags(Y(:,2), 0, I, I) + spdiags([X(2:I,2); 0], -1, I, I) + spdiags([0; Z(1:I-1,2)], 1, I, I);
        A = [A1, sparse(I,I); sparse(I,I), A2] + Aswitch;
    end

    B = (1/Delta + rho)*speye(2*I) - A;
    u_stacked = [u(:,1); u(:,2)];
    V_stacked = [V(:,1); V(:,2)];
    b = u_stacked + V_stacked/Delta;
    V_stacked = B\b;
    V = [V_stacked(1:I), V_stacked(I+1:2*I)];

    Vchange = V - v;
    v = V;
    if max(max(abs(Vchange))) < crit
        break
    end
end

% KFE
AT = A';
b_kfe = zeros(2*I, 1);
b_kfe(1) = 0.1;
AT(1,:) = [1, zeros(1, 2*I-1)];
gg = AT\b_kfe;
g_sum = gg'*ones(2*I,1)*da;
gg = gg/g_sum;
g = [gg(1:I), gg(I+1:2*I)];

% Aggregates
S = g(:,1)'*a*da + g(:,2)'*a*da;
L_F = da * (g(:,1)' * (z(1) * ell_F(:,1)) + g(:,2)' * (z(2) * ell_F(:,2)));
L_I = da * (g(:,1)' * (theta*z(1) * ell_I(:,1)) + g(:,2)' * (theta*z(2) * ell_I(:,2)));

v0_out = V;

end

