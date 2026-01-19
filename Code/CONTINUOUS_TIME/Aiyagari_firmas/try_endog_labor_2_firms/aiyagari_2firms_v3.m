%% AIYAGARI MODEL - 2 FIRMS (FORMAL/INFORMAL) - VERSION 3
% Based EXACTLY on codigo_prueba_v5.m with fzero for FOCs
% Uses lab_solve_2firms.m for coupled FOC system
% Author: Macroeconomic HA expert
% Date: 2026

clear all; clc; close all;

tic;

% =========================================================================
% 1. PARAMETERS
% =========================================================================

% Household
ga = 2;           % Risk aversion (γ)
rho = 0.05;       % Discount rate (ρ)
Frisch = 0.5;     % Frisch elasticity (φ)

% Two-Firm Parameters
psi_F = 1.0;      % Disutility formal
psi_I = 0.8;      % Disutility informal
theta = 0.5;      % Informal productivity attenuation

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
A_I = 0.15;       % w_I = A_I

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

maxit = 100;
crit = 10^(-6);
Delta = 1000;

dVf = zeros(I, 2);
dVb = zeros(I, 2);

Aswitch = [-speye(I)*la(1), speye(I)*la(1);
    speye(I)*la(2), -speye(I)*la(2)];

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
% 4. OUTER LOOP
% =========================================================================

fprintf('2-Firm Aiyagari v3 with fzero...\n');
fprintf('theta=%.2f, A_I=%.2f, psi_I=%.2f\n\n', theta, A_I, psi_I);

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

    % --- Precompute boundary labor using fzero (like Moll) ---
    ell_F_min = zeros(1,2);
    ell_I_min = zeros(1,2);
    c_min = zeros(1,2);
    dV_min = zeros(1,2);
    for j = 1:2
        params = [amin, z(j), w_F, w_I, theta, r, ga, Frisch, psi_F, psi_I];
        c_guess = w_F*z(j) + w_I*theta*z(j) + max(r*amin, 0.01);
        c_min(j) = fzero(@(c) lab_solve_2firms(c, params), c_guess, optimset('Display','off','TolFun',1e-10));
        c_min(j) = max(c_min(j), 1e-10);
        u_c = c_min(j)^(-ga);
        ell_F_min(j) = ((u_c * w_F * z(j)) / psi_F)^Frisch;
        ell_I_min(j) = ((u_c * w_I * theta * z(j)) / psi_I)^Frisch;
        dV_min(j) = c_min(j)^(-ga);
    end

    % =====================================================================
    % 5. INNER LOOP: HJB
    % =====================================================================

    for n = 1:maxit
        V = v;

        % --- 5.1. Finite Differences ---
        dVf(1:I-1,:) = (V(2:I,:) - V(1:I-1,:))/da;
        dVb(2:I,:) = (V(2:I,:) - V(1:I-1,:))/da;

        % --- 5.2. Boundary Conditions (using fzero like Moll) ---
        % At amax: solve FOC system
        for j = 1:2
            params = [amax, z(j), w_F, w_I, theta, r, ga, Frisch, psi_F, psi_I];
            c_guess = w_F*z(j) + w_I*theta*z(j) + max(r*amax, 0.01);
            c_up = fzero(@(c) lab_solve_2firms(c, params), c_guess, optimset('Display','off','TolFun',1e-10));
            c_up = max(c_up, 1e-10);
            u_c = c_up^(-ga);
            ell_F_up = ((u_c * w_F * z(j)) / psi_F)^Frisch;
            ell_I_up = ((u_c * w_I * theta * z(j)) / psi_I)^Frisch;
            income_up = w_F*z(j)*ell_F_up + w_I*theta*z(j)*ell_I_up + r*amax;
            dVf(I,j) = max(income_up, 1e-10)^(-ga);
        end

        % At amin: enforce state constraint
        dVb(1,:) = max((V(2,:) - V(1,:))/da, dV_min);

        % --- 5.3. UPWIND SCHEME WITH fzero (EXACT COPY of codigo_prueba_v5) ---

        % Forward difference
        dVf_pos = max(dVf, 1e-10);
        cf = real(dVf_pos .^ (-1/ga));
        ell_Ff = zeros(I, 2);
        ell_If = zeros(I, 2);
        for j = 1:2
            for i = 1:I
                % From FOC with known marginal utility u_c = dVf
                u_c = dVf_pos(i,j);
                ell_Ff(i,j) = ((u_c * w_F * z(j)) / psi_F)^Frisch;
                ell_If(i,j) = ((u_c * w_I * theta * z(j)) / psi_I)^Frisch;
            end
        end
        ssf = w_F*zz.*ell_Ff + w_I*theta*zz.*ell_If + r*aa - cf;

        % Backward difference
        dVb_pos = max(dVb, 1e-10);
        cb = real(dVb_pos .^ (-1/ga));
        ell_Fb = zeros(I, 2);
        ell_Ib = zeros(I, 2);
        for j = 1:2
            for i = 1:I
                u_c = dVb_pos(i,j);
                ell_Fb(i,j) = ((u_c * w_F * z(j)) / psi_F)^Frisch;
                ell_Ib(i,j) = ((u_c * w_I * theta * z(j)) / psi_I)^Frisch;
            end
        end
        ssb = w_F*zz.*ell_Fb + w_I*theta*zz.*ell_Ib + r*aa - cb;

        % Zero drift (using fzero through lab_solve_2firms)
        c0 = zeros(I, 2);
        ell_F0 = zeros(I, 2);
        ell_I0 = zeros(I, 2);
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

        % Upwind indicators
        If = ssf > 0;
        Ib = ssb < 0 & ~If;
        I0 = ~(If | Ib);

        % Upwind policies
        dV_Upwind = dVf_pos.*If + dVb_pos.*Ib + dV0.*I0;
        c = cf.*If + cb.*Ib + c0.*I0;
        ell_F = ell_Ff.*If + ell_Fb.*Ib + ell_F0.*I0;
        ell_I = ell_If.*If + ell_Ib.*Ib + ell_I0.*I0;

        % Utility with separate disutility
        u = c.^(1-ga)/(1-ga) ...
            - psi_F * ell_F.^(1+1/Frisch)/(1+1/Frisch) ...
            - psi_I * ell_I.^(1+1/Frisch)/(1+1/Frisch);

        % --- 5.4. Transition Matrix (EXACT COPY from codigo_prueba_v5) ---
        ss_Upwind = ssf .* If + ssb .* Ib + 0 .* I0;
        X = -min(ss_Upwind, 0) / da;
        Y = -max(ss_Upwind, 0) / da + min(ss_Upwind, 0) / da;
        Z = max(ss_Upwind, 0) / da;

        A1 = spdiags(Y(:,1), 0, I, I) + spdiags([X(2:I,1); 0], -1, I, I) + ...
            spdiags([0; Z(1:I-1,1)], 1, I, I);
        A2 = spdiags(Y(:,2), 0, I, I) + spdiags([X(2:I,2); 0], -1, I, I) + ...
            spdiags([0; Z(1:I-1,2)], 1, I, I);
        A = [A1, sparse(I,I); sparse(I,I), A2] + Aswitch;

        % Row sum adjustment
        row_sum = sum(A, 2);
        if max(abs(row_sum)) > 1e-9
            Y(:,1) = Y(:,1) - row_sum(1:I) / I;
            Y(:,2) = Y(:,2) - row_sum(I+1:2*I) / I;
            A1 = spdiags(Y(:,1), 0, I, I) + spdiags([X(2:I,1); 0], -1, I, I) + ...
                spdiags([0; Z(1:I-1,1)], 1, I, I);
            A2 = spdiags(Y(:,2), 0, I, I) + spdiags([X(2:I,2); 0], -1, I, I) + ...
                spdiags([0; Z(1:I-1,2)], 1, I, I);
            A = [A1, sparse(I,I); sparse(I,I), A2] + Aswitch;
            if max(abs(sum(A,2))) > 1e-9
                fprintf('Warning: Row sum failed at r=%.4f\n', r);
                break;
            end
        end

        % --- 5.5. Implicit solve ---
        B = (1/Delta + rho)*speye(2*I) - A;
        u_stacked = [u(:,1); u(:,2)];
        V_stacked = [V(:,1); V(:,2)];
        b = u_stacked + V_stacked/Delta;
        V_stacked = B\b;
        V = [V_stacked(1:I), V_stacked(I+1:2*I)];

        % Convergence
        Vchange = V - v;
        v = V;
        dist(n) = max(max(abs(Vchange)));

        if dist(n) < crit
            fprintf('Converged n=%d, r=%.4f\n', n, r);
            break
        end
    end

    % =====================================================================
    % 6. KFE
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

% Asset Supply and Demand
figure('Position', [100, 100, 800, 600])
Smax = max(S);
rrr = linspace(-0.06, 0.06, Ir);
KD_plot = (al*Aprod./(max(rrr + d, 1e-6))).^(1/(1 - al))*z_ave;

plot(S, r_grid, 'b-', 'LineWidth', 2.5)
hold on
plot(KD_plot, rrr, 'r-', 'LineWidth', 2.5)
plot([amin amin], [-0.06 0.06], 'k--')
plot([0 max(Smax,0.6)], [rho rho], 'k--')
plot([0 max(Smax,0.6)], [-d -d], 'k--')

xlabel('K', 'FontSize', 14)
ylabel('r', 'FontSize', 14)
title('Asset Supply and Demand - 2 Firms v3', 'FontSize', 16)
legend('S(r)', 'K^D(r)', 'Location', 'best')
grid on
xlim([0 max(Smax*1.1, 0.6)])
set(gca, 'FontSize', 12)
print('-dpng', [output_dir '/equilibrium_v3.png'], '-r300')

% Labor by Sector
figure('Position', [100, 100, 800, 600])
subplot(2,1,1)
plot(r_grid, L_F, 'b-', r_grid, L_I, 'r--', 'LineWidth', 2)
xlabel('r'); ylabel('Labor')
legend('L_F', 'L_I'); title('Aggregate Labor')
grid on

subplot(2,1,2)
plot(r_grid, L_I./(L_F + L_I + 1e-10), 'k-', 'LineWidth', 2)
xlabel('r'); ylabel('Share')
title('Informal Share')
grid on
print('-dpng', [output_dir '/labor_v3.png'], '-r300')

% Policy Functions
ir_mid = round(Ir/2);
figure('Position', [100, 100, 1000, 400])
subplot(1,2,1)
plot(a, ell_F_r(:,1,ir_mid), 'b-', a, ell_F_r(:,2,ir_mid), 'b--', 'LineWidth', 2)
xlabel('a'); ylabel('\ell_F')
legend('z_1', 'z_2'); title('Formal Labor')
grid on

subplot(1,2,2)
plot(a, ell_I_r(:,1,ir_mid), 'r-', a, ell_I_r(:,2,ir_mid), 'r--', 'LineWidth', 2)
xlabel('a'); ylabel('\ell_I')
legend('z_1', 'z_2'); title('Informal Labor')
grid on
print('-dpng', [output_dir '/policies_v3.png'], '-r300')

% =========================================================================
% 8. SAVE
% =========================================================================

save('aiyagari_2firms_v3.mat', ...
    'r_grid', 'S', 'KD', 'w_F_r', 'w_I_r', 'V_r', 'g_r', 'c_r', ...
    'ell_F_r', 'ell_I_r', 'L_F', 'L_I', 'a', 'z', ...
    'theta', 'A_I', 'psi_F', 'psi_I')

fprintf('\n=== COMPLETE ===\n')
