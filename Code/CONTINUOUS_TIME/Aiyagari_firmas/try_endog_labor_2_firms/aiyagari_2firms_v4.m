%% AIYAGARI MODEL - 2 FIRMS (FORMAL/INFORMAL) - VERSION 4
% With EQUILIBRIUM MODE selector:
%   EQUILIBRIUM_MODE = 1 → Partial Equilibrium (asset supply curve)
%   EQUILIBRIUM_MODE = 2 → General Equilibrium (find r* via bisection)
%
% Based on aiyagari_2firms_v3.m
% Author: Macroeconomic HA expert
% Date: 2026

clear all; clc; close all;

tic;

% =========================================================================
% 0. EQUILIBRIUM MODE SELECTOR
% =========================================================================

EQUILIBRIUM_MODE = 1;  % 1 = Partial (curve S(r)), 2 = General (find r*)

% =========================================================================
% 1. PARAMETERS
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

% =========================================================================
% 3. EQUILIBRIUM SOLVER
% =========================================================================

if EQUILIBRIUM_MODE == 1
    % =====================================================================
    % PARTIAL EQUILIBRIUM: Compute S(r) curve
    % =====================================================================

    Ir = 100;
    rmin = -0.0499;
    rmax = 0.049;
    r_grid = linspace(rmin, rmax, Ir);

    fprintf('=== PARTIAL EQUILIBRIUM MODE ===\n');
    fprintf('Computing asset supply curve S(r)...\n\n');

    % Initial guess for first r
    r = r_grid(1);
    KD_init = (al*Aprod/(r + d))^(1/(1 - al))*z_ave;
    w_F = (1 - al)*Aprod*(KD_init/z_ave)^al;
    w_I = A_I;
    v0(:,1) = (w_F*z(1) + w_I*theta*z(1) + max(r,0.01)*a).^(1-ga)/(1-ga)/rho;
    v0(:,2) = (w_F*z(2) + w_I*theta*z(2) + max(r,0.01)*a).^(1-ga)/(1-ga)/rho;

    for ir = 1:Ir
        r = r_grid(ir);
        [S(ir), KD(ir), w_F_r(ir), L_F(ir), L_I(ir), V, g, c, ell_F, ell_I, v0] = ...
            solve_partial_eq(r, v0, a, z, la, ga, rho, Frisch, psi_F, psi_I, theta, ...
            Aprod, al, d, A_I, z_ave, I, da, aa, zz, maxit, crit, Delta, Aswitch);

        % Store for later
        V_r(:,:,ir) = V;
        g_r(:,:,ir) = g;
        c_r(:,:,ir) = c;
        ell_F_r(:,:,ir) = ell_F;
        ell_I_r(:,:,ir) = ell_I;
        w_I_r(ir) = A_I;

        if mod(ir, 10) == 0
            fprintf('ir=%d/%d, S=%.4f, L_F=%.4f, L_I=%.4f\n', ir, Ir, S(ir), L_F(ir), L_I(ir));
        end
    end

elseif EQUILIBRIUM_MODE == 2
    % =====================================================================
    % GENERAL EQUILIBRIUM: Find r* via bisection
    % =====================================================================

    fprintf('=== GENERAL EQUILIBRIUM MODE ===\n');
    fprintf('Finding equilibrium r* via bisection...\n\n');

    % Bisection parameters
    r_low = -0.04;
    r_high = 0.045;
    tol_r = 1e-5;
    max_bisect = 50;

    % Initial guess
    r = (r_low + r_high)/2;
    KD_init = (al*Aprod/(r + d))^(1/(1 - al))*z_ave;
    w_F = (1 - al)*Aprod*(KD_init/z_ave)^al;
    w_I = A_I;
    v0(:,1) = (w_F*z(1) + w_I*theta*z(1) + max(r,0.01)*a).^(1-ga)/(1-ga)/rho;
    v0(:,2) = (w_F*z(2) + w_I*theta*z(2) + max(r,0.01)*a).^(1-ga)/(1-ga)/rho;

    % Bisection loop
    for iter = 1:max_bisect
        r = (r_low + r_high)/2;

        [S_mid, KD_mid, w_F_mid, L_F_mid, L_I_mid, V, g, c, ell_F, ell_I, v0] = ...
            solve_partial_eq(r, v0, a, z, la, ga, rho, Frisch, psi_F, psi_I, theta, ...
            Aprod, al, d, A_I, z_ave, I, da, aa, zz, maxit, crit, Delta, Aswitch);

        excess = S_mid - KD_mid;  % Excess supply

        fprintf('Bisect iter %2d: r=%.6f, S=%.4f, K^D=%.4f, excess=%.4f\n', ...
            iter, r, S_mid, KD_mid, excess);

        if abs(excess) < tol_r || (r_high - r_low) < tol_r
            fprintf('\n*** EQUILIBRIUM FOUND ***\n');
            break;
        end

        if excess > 0
            r_high = r;  % Too much supply → lower r
        else
            r_low = r;   % Too much demand → raise r
        end
    end

    % Store equilibrium values
    r_star = r;
    S_star = S_mid;
    K_star = KD_mid;
    w_F_star = w_F_mid;
    w_I_star = A_I;
    L_F_star = L_F_mid;
    L_I_star = L_I_mid;

    % For compatibility with plots
    Ir = 1;
    r_grid = r_star;
    S = S_star;
    KD = K_star;
    L_F = L_F_star;
    L_I = L_I_star;
    V_r(:,:,1) = V;
    g_r(:,:,1) = g;
    c_r(:,:,1) = c;
    ell_F_r(:,:,1) = ell_F;
    ell_I_r(:,:,1) = ell_I;
    w_F_r = w_F_star;
    w_I_r = w_I_star;
end

toc;

% =========================================================================
% 4. DISPLAY RESULTS
% =========================================================================

fprintf('\n========================================\n');
if EQUILIBRIUM_MODE == 1
    fprintf('PARTIAL EQUILIBRIUM RESULTS\n');
    fprintf('S(r) computed for %d values of r\n', Ir);
    fprintf('r range: [%.4f, %.4f]\n', rmin, rmax);
else
    fprintf('GENERAL EQUILIBRIUM RESULTS\n');
    fprintf('Equilibrium interest rate: r* = %.6f\n', r_star);
    fprintf('Equilibrium capital:       K* = %.4f\n', K_star);
    fprintf('Formal wage:               w_F* = %.4f\n', w_F_star);
    fprintf('Informal wage:             w_I* = %.4f\n', w_I_star);
    fprintf('Formal labor:              L_F* = %.4f\n', L_F_star);
    fprintf('Informal labor:            L_I* = %.4f\n', L_I_star);
    fprintf('Informal share:            %.2f%%\n', 100*L_I_star/(L_F_star+L_I_star));
end
fprintf('========================================\n\n');

% =========================================================================
% 5. PLOTS (Moll-style, each saved separately)
% =========================================================================

output_dir = 'output_graphs';
if ~exist(output_dir, 'dir'), mkdir(output_dir); end

% Get wages at equilibrium
if EQUILIBRIUM_MODE == 2
    w_F = w_F_star;
    w_I = w_I_star;
    r_eq = r_star;
else
    % Use mid-point of grid for partial equilibrium plots
    ir_mid = round(Ir/2);
    r_eq = r_grid(ir_mid);
    KD_mid = (al*Aprod/(r_eq + d))^(1/(1 - al))*z_ave;
    w_F = (1 - al)*Aprod*(KD_mid/z_ave)^al;
    w_I = A_I;
    g = g_r(:,:,ir_mid);
    c = c_r(:,:,ir_mid);
    ell_F = ell_F_r(:,:,ir_mid);
    ell_I = ell_I_r(:,:,ir_mid);
    V = V_r(:,:,ir_mid);
end

% --- Compute saving policy: s(a,z) = income - consumption ---
% Total income = w_F*z*ell_F + w_I*theta*z*ell_I + r*a
income = w_F*zz.*ell_F + w_I*theta*zz.*ell_I + r_eq*aa;
adot = income - c;  % Saving = income - consumption

% Plot limits (like Moll)
amax_plot = 5;
amin_plot = amin - 0.1;

% -------------------------------------------------------------------------
% FIGURE 1: SAVING POLICY s(a,z) - Moll style
% -------------------------------------------------------------------------
figure('Position', [100, 100, 700, 500])
plot(a, adot(:,1), 'b-', 'LineWidth', 2.5)
hold on
plot(a, adot(:,2), 'r--', 'LineWidth', 2.5)
plot(linspace(amin_plot, amax_plot, 100), zeros(1,100), 'k--', 'LineWidth', 1)
line([amin amin], [min(min(adot))*1.1 max(max(adot))*1.1], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1)

legend('s_1(a) - Low Productivity', 's_2(a) - High Productivity', 'Location', 'NorthEast')
xlabel('Wealth, $a$', 'FontSize', 14, 'Interpreter', 'latex')
ylabel('Savings, $s_i(a)$', 'FontSize', 14, 'Interpreter', 'latex')
title(sprintf('Saving Policy at r=%.4f', r_eq), 'FontSize', 14)
xlim([amin_plot amax_plot])
% Use simple ylim based on visible range
idx_max = min(round(amax_plot/da)+1, I);
ylim([min(min(adot(1:idx_max,:)))*1.2 max(max(adot(1:idx_max,:)))*1.2])
grid on
set(gca, 'FontSize', 12)
print('-dpng', [output_dir '/savings_policy.png'], '-r300')

% -------------------------------------------------------------------------
% FIGURE 2: WEALTH DISTRIBUTION g(a,z) - Moll style
% -------------------------------------------------------------------------
figure('Position', [100, 100, 700, 500])
plot(a, g(:,1), 'b-', 'LineWidth', 2.5)
hold on
plot(a, g(:,2), 'r--', 'LineWidth', 2.5)
line([amin amin], [0 max(max(g))*1.1], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1)

legend('g_1(a) - Low Productivity', 'g_2(a) - High Productivity', 'Location', 'NorthEast')
xlabel('Wealth, $a$', 'FontSize', 14, 'Interpreter', 'latex')
ylabel('Density, $g_i(a)$', 'FontSize', 14, 'Interpreter', 'latex')
title(sprintf('Stationary Distribution at r=%.4f', r_eq), 'FontSize', 14)
xlim([amin_plot amax_plot])
ylim([0 max(max(g(1:idx_max,:)))*1.1])
grid on
set(gca, 'FontSize', 12)
print('-dpng', [output_dir '/wealth_distribution.png'], '-r300')

% -------------------------------------------------------------------------
% FIGURE 3: CONSUMPTION POLICY c(a,z)
% -------------------------------------------------------------------------
figure('Position', [100, 100, 700, 500])
plot(a, c(:,1), 'b-', 'LineWidth', 2.5)
hold on
plot(a, c(:,2), 'r--', 'LineWidth', 2.5)

legend('c_1(a) - Low Productivity', 'c_2(a) - High Productivity', 'Location', 'SouthEast')
xlabel('Wealth, $a$', 'FontSize', 14, 'Interpreter', 'latex')
ylabel('Consumption, $c_i(a)$', 'FontSize', 14, 'Interpreter', 'latex')
title(sprintf('Consumption Policy at r=%.4f', r_eq), 'FontSize', 14)
xlim([amin_plot amax_plot])
grid on
set(gca, 'FontSize', 12)
print('-dpng', [output_dir '/consumption_policy.png'], '-r300')

% -------------------------------------------------------------------------
% FIGURE 4: FORMAL LABOR POLICY ell_F(a,z)
% -------------------------------------------------------------------------
figure('Position', [100, 100, 700, 500])
plot(a, ell_F(:,1), 'b-', 'LineWidth', 2.5)
hold on
plot(a, ell_F(:,2), 'b--', 'LineWidth', 2.5)

legend('ell_F(z_1)', 'ell_F(z_2)', 'Location', 'NorthEast')
xlabel('Wealth, $a$', 'FontSize', 14, 'Interpreter', 'latex')
ylabel('Formal Labor, $\ell_F(a)$', 'FontSize', 14, 'Interpreter', 'latex')
title(sprintf('Formal Labor Policy at r=%.4f', r_eq), 'FontSize', 14)
xlim([amin_plot amax_plot])
grid on
set(gca, 'FontSize', 12)
print('-dpng', [output_dir '/labor_formal.png'], '-r300')

% -------------------------------------------------------------------------
% FIGURE 5: INFORMAL LABOR POLICY ell_I(a,z)
% -------------------------------------------------------------------------
figure('Position', [100, 100, 700, 500])
plot(a, ell_I(:,1), 'r-', 'LineWidth', 2.5)
hold on
plot(a, ell_I(:,2), 'r--', 'LineWidth', 2.5)

legend('ell_I(z_1)', 'ell_I(z_2)', 'Location', 'NorthEast')
xlabel('Wealth, $a$', 'FontSize', 14, 'Interpreter', 'latex')
ylabel('Informal Labor, $\ell_I(a)$', 'FontSize', 14, 'Interpreter', 'latex')
title(sprintf('Informal Labor Policy at r=%.4f', r_eq), 'FontSize', 14)
xlim([amin_plot amax_plot])
grid on
set(gca, 'FontSize', 12)
print('-dpng', [output_dir '/labor_informal.png'], '-r300')

% -------------------------------------------------------------------------
% FIGURE 6: VALUE FUNCTION V(a,z)
% -------------------------------------------------------------------------
figure('Position', [100, 100, 700, 500])
plot(a, V(:,1), 'b-', 'LineWidth', 2.5)
hold on
plot(a, V(:,2), 'r--', 'LineWidth', 2.5)

legend('V_1(a) - Low Productivity', 'V_2(a) - High Productivity', 'Location', 'SouthEast')
xlabel('Wealth, $a$', 'FontSize', 14, 'Interpreter', 'latex')
ylabel('Value Function, $V_i(a)$', 'FontSize', 14, 'Interpreter', 'latex')
title(sprintf('Value Function at r=%.4f', r_eq), 'FontSize', 14)
xlim([amin_plot amax_plot])
grid on
set(gca, 'FontSize', 12)
print('-dpng', [output_dir '/value_function.png'], '-r300')

% -------------------------------------------------------------------------
% FIGURE 7: EQUILIBRIUM (Partial or General)
% -------------------------------------------------------------------------
if EQUILIBRIUM_MODE == 1
    % Partial: plot S(r) vs K^D(r)
    figure('Position', [100, 100, 700, 500])
    Smax = max(S);
    plot(S, r_grid, 'b-', 'LineWidth', 2.5)
    hold on
    % Use the computed KD (which accounts for endogenous labor) instead of synthetic one
    plot(KD, r_grid, 'r-', 'LineWidth', 2.5)
    plot([0 max(Smax,0.6)], [rho rho], 'k--', 'LineWidth', 1)
    plot([0 max(Smax,0.6)], [-d -d], 'k--', 'LineWidth', 1)

    xlabel('Capital, $K$', 'FontSize', 14, 'Interpreter', 'latex')
    ylabel('Interest Rate, $r$', 'FontSize', 14, 'Interpreter', 'latex')
    title('Partial Equilibrium: Asset Supply and Demand', 'FontSize', 14)
    legend('S(r) - Supply', 'K^D(r) - Demand', 'Location', 'best')
    xlim([0 max(Smax*1.1, 0.6)])
    grid on
    set(gca, 'FontSize', 12)
    print('-dpng', [output_dir '/equilibrium_partial.png'], '-r300')
else
    % General: show equilibrium point
    figure('Position', [100, 100, 700, 500])
    rrr = linspace(-0.06, 0.06, 100);
    % Update: Plot KD curve assuming Labor stays at equilibrium level L_F_star
    % This shows the demand curve CONSISTENT with the final labor supply.
    KD_plot = (al*Aprod./(max(rrr + d, 1e-6))).^(1/(1 - al)) * L_F_star;

    plot(KD_plot, rrr, 'r-', 'LineWidth', 2.5)
    hold on
    plot(K_star, r_star, 'go', 'MarkerSize', 15, 'MarkerFaceColor', 'g', 'LineWidth', 2)
    plot([0 K_star*2], [r_star r_star], 'k--', 'LineWidth', 1)
    plot([K_star K_star], [-0.06 r_star], 'k--', 'LineWidth', 1)

    xlabel('Capital, $K$', 'FontSize', 14, 'Interpreter', 'latex')
    ylabel('Interest Rate, $r$', 'FontSize', 14, 'Interpreter', 'latex')
    title(sprintf('General Equilibrium: $r^*=%.4f$, $K^*=%.3f$', r_star, K_star), 'FontSize', 14, 'Interpreter', 'latex')
    legend('K^D(r)', sprintf('Equilibrium (r*=%.3f%%)', r_star*100), 'Location', 'NorthEast')
    xlim([0 K_star*2])
    ylim([-0.06 0.06])
    grid on
    set(gca, 'FontSize', 12)
    print('-dpng', [output_dir '/equilibrium_general.png'], '-r300')
end

% -------------------------------------------------------------------------
% FIGURE 8: INFORMAL LABOR SHARE by wealth
% -------------------------------------------------------------------------
figure('Position', [100, 100, 700, 500])
informal_share = ell_I ./ (ell_F + ell_I + 1e-10);
plot(a, informal_share(:,1), 'k-', 'LineWidth', 2.5)
hold on
plot(a, informal_share(:,2), 'k--', 'LineWidth', 2.5)

legend('Low Prod.', 'High Prod.', 'Location', 'best')
xlabel('Wealth, $a$', 'FontSize', 14, 'Interpreter', 'latex')
ylabel('Informal Share, $\ell_I / (\ell_F + \ell_I)$', 'FontSize', 14, 'Interpreter', 'latex')
title('Informal Labor Share by Wealth', 'FontSize', 14)
xlim([amin_plot amax_plot])
ylim([0 1])
grid on
set(gca, 'FontSize', 12)
print('-dpng', [output_dir '/informal_share.png'], '-r300')

% -------------------------------------------------------------------------
% FIGURE 9: INCOME DECOMPOSITION (Formal vs Informal)
% -------------------------------------------------------------------------
% Since there's one consumption c(a,z), we show income by source instead
income_formal = w_F * zz .* ell_F;      % Formal labor income
income_informal = w_I * theta * zz .* ell_I;  % Informal labor income
income_assets = r_eq * aa;               % Asset income

figure('Position', [100, 100, 800, 500])
subplot(1,2,1)
plot(a, income_formal(:,1), 'b-', 'LineWidth', 2)
hold on
plot(a, income_informal(:,1), 'r--', 'LineWidth', 2)
plot(a, income_assets(:,1), 'k:', 'LineWidth', 2)
legend('Formal (w_F z ell_F)', 'Informal (w_I theta z ell_I)', 'Assets (r a)', 'Location', 'best')
xlabel('Wealth, $a$', 'FontSize', 12, 'Interpreter', 'latex')
ylabel('Income', 'FontSize', 12)
title('Income by Source (Low Productivity)', 'FontSize', 14)
xlim([amin_plot amax_plot])
grid on

subplot(1,2,2)
plot(a, income_formal(:,2), 'b-', 'LineWidth', 2)
hold on
plot(a, income_informal(:,2), 'r--', 'LineWidth', 2)
plot(a, income_assets(:,2), 'k:', 'LineWidth', 2)
legend('Formal', 'Informal', 'Assets', 'Location', 'best')
xlabel('Wealth, $a$', 'FontSize', 12, 'Interpreter', 'latex')
ylabel('Income', 'FontSize', 12)
title('Income by Source (High Productivity)', 'FontSize', 14)
xlim([amin_plot amax_plot])
grid on
print('-dpng', [output_dir '/income_decomposition.png'], '-r300')

% -------------------------------------------------------------------------
% FIGURE 10: FORMAL vs INFORMAL LABOR INCOME SHARE
% -------------------------------------------------------------------------
total_labor_income = income_formal + income_informal;
formal_income_share = income_formal ./ (total_labor_income + 1e-10);

figure('Position', [100, 100, 700, 500])
plot(a, formal_income_share(:,1), 'b-', 'LineWidth', 2.5)
hold on
plot(a, formal_income_share(:,2), 'b--', 'LineWidth', 2.5)

legend('Low Prod.', 'High Prod.', 'Location', 'best')
xlabel('Wealth, $a$', 'FontSize', 14, 'Interpreter', 'latex')
ylabel('Formal Income Share', 'FontSize', 14)
title('Formal Labor Income / Total Labor Income', 'FontSize', 14)
xlim([amin_plot amax_plot])
ylim([0 1])
grid on
set(gca, 'FontSize', 12)
print('-dpng', [output_dir '/formal_income_share.png'], '-r300')

% -------------------------------------------------------------------------
% FIGURE 11: AGGREGATE INCOME DISTRIBUTION (Bar chart with std)
% -------------------------------------------------------------------------
% Calculate aggregate income by source using g(a,z) as weights

% Aggregate over both productivity states
g_total = g(:,1) + g(:,2);  % Combined density
g_total = g_total / (sum(g_total)*da);  % Normalize

% Mean income by source (weighted by g)
mean_formal = da * (g(:,1)' * income_formal(:,1) + g(:,2)' * income_formal(:,2));
mean_informal = da * (g(:,1)' * income_informal(:,1) + g(:,2)' * income_informal(:,2));
mean_assets = da * (g(:,1)' * income_assets(:,1) + g(:,2)' * income_assets(:,2));
mean_total = mean_formal + mean_informal + mean_assets;

% Variance calculation
var_formal = da * (g(:,1)' * (income_formal(:,1) - mean_formal).^2 + g(:,2)' * (income_formal(:,2) - mean_formal).^2);
var_informal = da * (g(:,1)' * (income_informal(:,1) - mean_informal).^2 + g(:,2)' * (income_informal(:,2) - mean_informal).^2);
var_assets = da * (g(:,1)' * (income_assets(:,1) - mean_assets).^2 + g(:,2)' * (income_assets(:,2) - mean_assets).^2);
total_income = income_formal + income_informal + income_assets;
mean_total_income = da * (g(:,1)' * total_income(:,1) + g(:,2)' * total_income(:,2));
var_total = da * (g(:,1)' * (total_income(:,1) - mean_total_income).^2 + g(:,2)' * (total_income(:,2) - mean_total_income).^2);

std_formal = sqrt(var_formal);
std_informal = sqrt(var_informal);
std_assets = sqrt(var_assets);
std_total = sqrt(var_total);

% Bar chart
figure('Position', [100, 100, 800, 500])
means = [mean_formal, mean_informal, mean_assets, mean_total_income];
stds = [std_formal, std_informal, std_assets, std_total];

bar_handle = bar(1:4, means, 0.6, 'FaceColor', 'flat');
bar_handle.CData(1,:) = [0 0.4470 0.7410];  % Blue for formal
bar_handle.CData(2,:) = [0.8500 0.3250 0.0980];  % Red for informal
bar_handle.CData(3,:) = [0.5 0.5 0.5];  % Gray for assets
bar_handle.CData(4,:) = [0.4660 0.6740 0.1880];  % Green for total

hold on
errorbar(1:4, means, stds, 'k.', 'LineWidth', 2, 'CapSize', 10)

set(gca, 'XTickLabel', {'Formal', 'Informal', 'Assets', 'Total'}, 'FontSize', 12)
ylabel('Income', 'FontSize', 14)
title('Aggregate Income Distribution by Source', 'FontSize', 14)
grid on

% Add value labels on bars
for i = 1:4
    text(i, means(i) + stds(i) + 0.01, sprintf('%.3f', means(i)), ...
        'HorizontalAlignment', 'center', 'FontSize', 10)
end

% Add legend explaining
text(0.5, max(means)*1.3, sprintf('Mean \\pm Std'), 'FontSize', 12, 'FontWeight', 'bold')
print('-dpng', [output_dir '/income_distribution_bars.png'], '-r300')

% -------------------------------------------------------------------------
% FIGURE 12: INCOME SHARE PIE CHART
% -------------------------------------------------------------------------
figure('Position', [100, 100, 600, 500])
shares = [mean_formal, mean_informal, mean_assets];
labels = {sprintf('Formal (%.1f%%)', 100*mean_formal/mean_total_income), ...
    sprintf('Informal (%.1f%%)', 100*mean_informal/mean_total_income), ...
    sprintf('Assets (%.1f%%)', 100*mean_assets/mean_total_income)};
pie(shares, labels)
title('Income Composition at Equilibrium', 'FontSize', 14)
colormap([0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.5 0.5 0.5])
print('-dpng', [output_dir '/income_pie_chart.png'], '-r300')

fprintf('Graphs saved to %s/\n', output_dir);

% =========================================================================
% 6. SAVE
% =========================================================================

if EQUILIBRIUM_MODE == 1
    save('aiyagari_2firms_partial.mat', ...
        'r_grid', 'S', 'KD', 'w_F_r', 'w_I_r', 'V_r', 'g_r', 'c_r', ...
        'ell_F_r', 'ell_I_r', 'L_F', 'L_I', 'a', 'z', ...
        'theta', 'A_I', 'psi_F', 'psi_I', 'EQUILIBRIUM_MODE')
else
    save('aiyagari_2firms_general.mat', ...
        'r_star', 'K_star', 'w_F_star', 'w_I_star', 'L_F_star', 'L_I_star', ...
        'V', 'g', 'c', 'ell_F', 'ell_I', 'a', 'z', ...
        'theta', 'A_I', 'psi_F', 'psi_I', 'EQUILIBRIUM_MODE')
end

fprintf('=== COMPLETE ===\n')


% =========================================================================
% HELPER FUNCTION: Solve partial equilibrium for given r
% =========================================================================

function [S, KD, w_F, L_F, L_I, V, g, c, ell_F, ell_I, v0_out] = ...
    solve_partial_eq(r, v0, a, z, la, ga, rho, Frisch, psi_F, psi_I, theta, ...
    Aprod, al, d, A_I, z_ave, I, da, aa, zz, maxit, crit, Delta, Aswitch)

% Firm prices
% 1. Capital-Labor ratio determined by r (FOC)
k_ratio = (al*Aprod/(r + d))^(1/(1 - al));
% 2. Wage determined by k_ratio
w_F = (1 - al)*Aprod*(k_ratio)^al;
w_I = A_I;

v = v0;
dVf = zeros(I, 2);
dVb = zeros(I, 2);
amin = a(1);
amax = a(end);

% Precompute boundary
dV_min = zeros(1,2);
for j = 1:2
    params = [amin, z(j), w_F, w_I, theta, r, ga, Frisch, psi_F, psi_I];
    c_guess = w_F*z(j) + w_I*theta*z(j) + max(r*amin, 0.01);
    c_min = fzero(@(c) lab_solve_2firms(c, params), c_guess, optimset('Display','off','TolFun',1e-10));
    dV_min(j) = max(c_min, 1e-10)^(-ga);
end

% HJB iteration
for n = 1:maxit
    V = v;

    dVf(1:I-1,:) = (V(2:I,:) - V(1:I-1,:))/da;
    dVb(2:I,:) = (V(2:I,:) - V(1:I-1,:))/da;

    % Boundaries
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

    % Zero drift
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

% UPDATE KD consistent with Endogenous Labor
% The firm FOC implies K/L ratio is fixed by r.
% So Capital Demand must scale with Labor Supply L_F.
KD = k_ratio * L_F;

v0_out = V;

end
