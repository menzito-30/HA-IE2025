%% AIYAGARI MODEL - 2 FIRMS (FORMAL/INFORMAL) - VERSION 2
% Based on codigo_prueba_v5.m (1 firma) - MINIMAL CHANGES for 2 firms
% Author: Adapted for formal/informal duality
% Date: 2026

clear all; clc; close all;

tic;

% =========================================================================
% 1. PARAMETERS (Same as 1-firma code + 2-firm additions)
% =========================================================================

% --- Household Parameters ---
ga = 2;           % Relative risk aversion (γ)
rho = 0.05;       % Subjective discount rate (ρ)
Frisch = 0.5;     % Frisch elasticity of labor supply (φ)

% --- TWO-FIRM PARAMETERS ---
psi_F = 1.0;      % Disutility weight formal (normalized to 1)
psi_I = 0.8;      % Disutility weight informal (< 1: easier work)
theta = 0.5;      % Informal sector productivity attenuation

% --- Productivity States (z) ---
z1 = 0.2;
z2 = 2*z1;
z = [z1, z2];
la1 = 1;
la2 = 1;
la = [la1, la2];
z_ave = (z1*la2 + z2*la1)/(la1 + la2);

% --- Formal Firm Parameters ---
Aprod = 0.3;
al = 1/3;
d = 0.05;

% --- Informal Firm Parameters ---
A_I = 0.15;       % w_I = A_I (constant)

% =========================================================================
% 2. GRIDS
% =========================================================================

I = 500;
amin = 0;
amax = 20;
a = linspace(amin, amax, I)';
da = (amax - amin)/(I - 1);

aa = [a, a];
zz = ones(I,1)*z;

% HJB iteration parameters
maxit = 100;
crit = 10^(-6);
Delta = 1000;

% Initialize
dVf = zeros(I, 2);
dVb = zeros(I, 2);

% Productivity transition matrix
Aswitch = [-speye(I)*la(1), speye(I)*la(1);
    speye(I)*la(2), -speye(I)*la(2)];

% Interest rate grid
Ir = 100;
rmin = -0.0499;
rmax = 0.049;
r_grid = linspace(rmin, rmax, Ir);

% =========================================================================
% 3. INITIAL GUESS
% =========================================================================
r = r_grid(1);
KD = (al*Aprod/(r + d))^(1/(1 - al))*z_ave;
w_F = (1 - al)*Aprod*(KD/z_ave)^al;
w_I = A_I;

v0(:,1) = (w_F*z(1) + w_I*theta*z(1) + max(r,0.01)*a).^(1-ga)/(1-ga)/rho;
v0(:,2) = (w_F*z(2) + w_I*theta*z(2) + max(r,0.01)*a).^(1-ga)/(1-ga)/rho;

% =========================================================================
% 4. OUTER LOOP: INTEREST RATES
% =========================================================================

fprintf('Starting 2-Firm Aiyagari Model v2...\n');
fprintf('Parameters: theta=%.2f, A_I=%.2f, psi_I=%.2f\n', theta, A_I, psi_I);

for ir = 1:Ir

    r = r_grid(ir);

    % Firm prices
    KD(ir) = (al*Aprod/(r + d))^(1/(1 - al))*z_ave;
    w_F = (1 - al)*Aprod*(KD(ir)/z_ave)^al;
    w_I = A_I;
    w_F_r(ir) = w_F;
    w_I_r(ir) = w_I;

    % Warm start
    if ir > 1
        v0 = V_r(:,:,ir - 1);
    end
    v = v0;

    % =====================================================================
    % 5. INNER LOOP: SOLVE HJB
    % =====================================================================

    for n = 1:maxit
        V = v;

        % --- 5.1. Finite Difference Approximations ---
        dVf(1:I-1,:) = (V(2:I,:) - V(1:I-1,:))/da;
        dVb(2:I,:) = (V(2:I,:) - V(1:I-1,:))/da;

        % Boundary at amax
        dVf(I,:) = (w_F*z + w_I*theta*z + r*amax).^(-ga);

        % Boundary at amin
        dVb(1,:) = (w_F*z + w_I*theta*z + r*amin).^(-ga);

        % Ensure positive (avoid complex numbers)
        dVf = max(dVf, 1e-10);
        dVb = max(dVb, 1e-10);

        % --- 5.2. POLICIES FROM FOC (ANALYTICAL SOLUTION) ---
        % From FOCs with SEPARABLE disutility:
        %   psi_F * ell_F^(1/phi) = u_c * w_F * z  =>  ell_F = ((u_c * w_F * z)/psi_F)^phi
        %   psi_I * ell_I^(1/phi) = u_c * w_I*theta*z  =>  ell_I = ((u_c * w_I*theta*z)/psi_I)^phi
        % where u_c = c^(-ga) = V'

        % Forward difference policies
        cf = dVf.^(-1/ga);
        % Labor from FOC (with c = dVf^(-1/ga), so dVf = c^(-ga) = u_c)
        ell_Ff = ((dVf .* w_F .* zz) / psi_F).^Frisch;
        ell_If = ((dVf .* w_I .* theta .* zz) / psi_I).^Frisch;
        ssf = w_F.*zz.*ell_Ff + w_I.*theta.*zz.*ell_If + r.*aa - cf;

        % Backward difference policies
        cb = dVb.^(-1/ga);
        ell_Fb = ((dVb .* w_F .* zz) / psi_F).^Frisch;
        ell_Ib = ((dVb .* w_I .* theta .* zz) / psi_I).^Frisch;
        ssb = w_F.*zz.*ell_Fb + w_I.*theta.*zz.*ell_Ib + r.*aa - cb;

        % Zero drift (c = income)
        % For zero drift, need to solve: c = w_F*z*ell_F + w_I*theta*z*ell_I + r*a
        % with FOCs. This gives a fixed point.
        % Approximation: use forward labor with zero-drift consumption
        c0 = w_F.*zz.*ell_Ff + w_I.*theta.*zz.*ell_If + r.*aa;
        c0 = max(c0, 1e-10);
        ell_F0 = ell_Ff;
        ell_I0 = ell_If;
        dV0 = c0.^(-ga);

        % --- 5.3. Upwind Scheme ---
        If = ssf > 0;
        Ib = ssb < 0 & ~If;
        I0 = ~(If | Ib);

        % Policies
        c = cf.*If + cb.*Ib + c0.*I0;
        ell_F = ell_Ff.*If + ell_Fb.*Ib + ell_F0.*I0;
        ell_I = ell_If.*If + ell_Ib.*Ib + ell_I0.*I0;

        % Ensure positive
        c = max(c, 1e-10);
        ell_F = max(ell_F, 0);
        ell_I = max(ell_I, 0);

        % Utility with SEPARATE disutility
        u = c.^(1-ga)/(1-ga) ...
            - psi_F * ell_F.^(1+1/Frisch)/(1+1/Frisch) ...
            - psi_I * ell_I.^(1+1/Frisch)/(1+1/Frisch);

        % --- 5.4. Transition Matrix A ---
        % EXACT COPY from codigo_prueba_v5.m
        ss_Upwind = ssf .* If + ssb .* Ib + 0 .* I0;
        X = -min(ss_Upwind, 0) / da;
        Y = -max(ss_Upwind, 0) / da + min(ss_Upwind, 0) / da;
        Z = max(ss_Upwind, 0) / da;

        A1 = spdiags(Y(:,1), 0, I, I) + spdiags([X(2:I,1); 0], -1, I, I) + ...
            spdiags([0; Z(1:I-1,1)], 1, I, I);
        A2 = spdiags(Y(:,2), 0, I, I) + spdiags([X(2:I,2); 0], -1, I, I) + ...
            spdiags([0; Z(1:I-1,2)], 1, I, I);
        A = [A1, sparse(I,I); sparse(I,I), A2] + Aswitch;

        % Check and adjust row sums (from codigo_prueba_v5.m)
        row_sum = sum(A, 2);
        if max(abs(row_sum)) > 1e-9
            Y(:,1) = Y(:,1) - row_sum(1:I) / I;
            Y(:,2) = Y(:,2) - row_sum(I+1:2*I) / I;
            A1 = spdiags(Y(:,1), 0, I, I) + spdiags([X(2:I,1); 0], -1, I, I) + ...
                spdiags([0; Z(1:I-1,1)], 1, I, I);
            A2 = spdiags(Y(:,2), 0, I, I) + spdiags([X(2:I,2); 0], -1, I, I) + ...
                spdiags([0; Z(1:I-1,2)], 1, I, I);
            A = [A1, sparse(I,I); sparse(I,I), A2] + Aswitch;
        end

        % --- 5.5. Solve Implicit System ---
        B = (1/Delta + rho)*speye(2*I) - A;
        u_stacked = [u(:,1); u(:,2)];
        V_stacked = [V(:,1); V(:,2)];
        b = u_stacked + V_stacked/Delta;
        V_stacked = B\b;
        V = [V_stacked(1:I), V_stacked(I+1:2*I)];

        % Convergence check
        Vchange = V - v;
        v = V;
        dist(n) = max(max(abs(Vchange)));

        if dist(n) < crit
            fprintf('Converged at n=%d for r=%.4f\n', n, r);
            break
        end
    end

    % =====================================================================
    % 6. KFE (EXACT COPY from codigo_prueba_v5.m)
    % =====================================================================

    AT = A';
    b = zeros(2*I, 1);
    i_fix = 1;
    b(i_fix) = 0.1;
    row = [zeros(1, i_fix-1), 1, zeros(1, 2*I - i_fix)];
    AT(i_fix,:) = row;

    gg = AT\b;
    g_sum = gg'*ones(2*I,1)*da;
    gg = gg/g_sum;

    g = [gg(1:I), gg(I+1:2*I)];

    % Store
    g_r(:,:,ir) = g;
    V_r(:,:,ir) = V;
    c_r(:,:,ir) = c;
    ell_F_r(:,:,ir) = ell_F;
    ell_I_r(:,:,ir) = ell_I;

    % Aggregates
    S(ir) = g(:,1)'*a*da + g(:,2)'*a*da;
    L_F(ir) = da * (g(:,1)' * (z(1) * ell_F(:,1)) + g(:,2)' * (z(2) * ell_F(:,2)));
    L_I(ir) = da * (g(:,1)' * (theta*z(1) * ell_I(:,1)) + g(:,2)' * (theta*z(2) * ell_I(:,2)));

    if mod(ir, 10) == 0
        fprintf('ir=%d/%d, S=%.4f, L_F=%.4f, L_I=%.4f\n', ir, Ir, S(ir), L_F(ir), L_I(ir));
    end
end

toc;

% =========================================================================
% 7. PLOTS
% =========================================================================

output_dir = 'output_graphs';
if ~exist(output_dir, 'dir'), mkdir(output_dir); end

% Figure 1: Asset Supply and Demand
figure('Position', [100, 100, 800, 600])
Smax = max(S);
amin1 = amin - 0.02;
aaa = linspace(amin1, max(Smax,0.6), Ir);
rrr = linspace(-0.06, 0.06, Ir);
KD_plot = (al*Aprod./(max(rrr + d, 1e-6))).^(1/(1 - al))*z_ave;

plot(S, r_grid, 'b-', 'LineWidth', 2.5)
hold on
plot(KD_plot, rrr, 'r-', 'LineWidth', 2.5)
plot([amin amin], [-0.06 0.06], 'k--', 'LineWidth', 1)
plot([amin1 max(Smax,0.6)], [rho rho], 'k--', 'LineWidth', 1)
plot([amin1 max(Smax,0.6)], [-d -d], 'k--', 'LineWidth', 1)

xlabel('$K$', 'FontSize', 16, 'interpreter', 'latex')
ylabel('$r$', 'FontSize', 16, 'interpreter', 'latex')
title('Asset Supply and Demand - 2 Firms v2', 'FontSize', 18)
legend('S(r)', 'K^D(r)', 'Location', 'best')
grid on
xlim([amin1 max(Smax*1.1, 0.6)])
set(gca, 'FontSize', 14)
print('-dpng', [output_dir '/equilibrium_2firms_v2.png'], '-r300')

% Figure 2: Labor by Sector
figure('Position', [100, 100, 800, 600])
subplot(2,1,1)
plot(r_grid, L_F, 'b-', 'LineWidth', 2)
hold on
plot(r_grid, L_I, 'r--', 'LineWidth', 2)
xlabel('r'); ylabel('Labor')
title('Aggregate Labor by Sector')
legend('L_F', 'L_I')
grid on

subplot(2,1,2)
plot(r_grid, L_I./(L_F + L_I + 1e-10), 'k-', 'LineWidth', 2)
xlabel('r'); ylabel('Informal Share')
title('L_I / (L_F + L_I)')
grid on
print('-dpng', [output_dir '/labor_by_sector_v2.png'], '-r300')

% Figure 3: Policy Functions
ir_mid = round(Ir/2);
figure('Position', [100, 100, 1000, 400])
subplot(1,2,1)
plot(a, ell_F_r(:,1,ir_mid), 'b-', a, ell_F_r(:,2,ir_mid), 'b--', 'LineWidth', 2)
xlabel('a'); ylabel('\ell_F')
legend('z_1', 'z_2')
title('Formal Labor')
grid on

subplot(1,2,2)
plot(a, ell_I_r(:,1,ir_mid), 'r-', a, ell_I_r(:,2,ir_mid), 'r--', 'LineWidth', 2)
xlabel('a'); ylabel('\ell_I')
legend('z_1', 'z_2')
title('Informal Labor')
grid on
print('-dpng', [output_dir '/labor_policies_v2.png'], '-r300')

% =========================================================================
% 8. SAVE
% =========================================================================

save('aiyagari_2firms_v2_results.mat', ...
    'r_grid', 'S', 'KD', 'w_F_r', 'w_I_r', 'V_r', 'g_r', 'c_r', ...
    'ell_F_r', 'ell_I_r', 'L_F', 'L_I', 'a', 'z', ...
    'theta', 'A_I', 'psi_F', 'psi_I')

fprintf('\n=== COMPLETE ===\n')
