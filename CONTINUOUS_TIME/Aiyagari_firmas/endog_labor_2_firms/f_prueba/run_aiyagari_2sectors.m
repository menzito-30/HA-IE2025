%% Main Script: Aiyagari Model with Two Sectors and Endogenous Labor Supply
% 
% This script solves the Aiyagari model with:
%   - Formal and informal sectors
%   - Endogenous labor supply (intensive margin)
%   - Proportional taxation on formal sector
%   - Balanced government budget
%
% Based on continuous time methods (Pontus Rendahl, Benjamin Moll)

clear; clc; close all;

%% Set parameters
params = struct();

% Preferences
params.rho = 0.05;           % Discount rate
params.gamma = 2.0;          % Risk aversion (CRRA)
params.frisch = 0.5;         % Frisch elasticity of labor supply

% Technology
params.alpha = 0.33;         % Capital share in production
params.delta = 0.05;         % Depreciation rate

% Productivity process (2-state Poisson)
params.z1 = 0.8;             % Low productivity
params.z2 = 1.2;             % High productivity
params.lambda1 = 0.5;        % Intensity z1 -> z2
params.lambda2 = 0.5;        % Intensity z2 -> z1

% Asset grid
params.abar = 0.0;           % Borrowing constraint (natural borrowing limit)
params.amax = 50;            % Maximum assets
params.N = 200;              % Grid points

% Government
params.tau_target = 0.20;    % Proportional tax on formal income (20%)

% Sectors
params.wI_wF_ratio = 0.75;   % Informal wage as fraction of formal wage

% Numerical parameters
params.Gamma = 1e10;         % Implicit method parameter (large = pure implicit)

%% Display parameters
fprintf('=================================================\n');
fprintf('   Aiyagari Model - Two Sectors with Labor Supply\n');
fprintf('=================================================\n\n');

fprintf('PREFERENCES:\n');
fprintf('  Discount rate (rho):       %.4f\n', params.rho);
fprintf('  Risk aversion (gamma):     %.4f\n', params.gamma);
fprintf('  Frisch elasticity:         %.4f\n\n', params.frisch);

fprintf('TECHNOLOGY:\n');
fprintf('  Capital share (alpha):     %.4f\n', params.alpha);
fprintf('  Depreciation (delta):      %.4f\n\n', params.delta);

fprintf('PRODUCTIVITY:\n');
fprintf('  Low state (z1):            %.4f\n', params.z1);
fprintf('  High state (z2):           %.4f\n', params.z2);
fprintf('  Transition rate (lambda1): %.4f\n', params.lambda1);
fprintf('  Transition rate (lambda2): %.4f\n\n', params.lambda2);

fprintf('GOVERNMENT:\n');
fprintf('  Tax rate (tau):            %.4f\n\n', params.tau_target);

fprintf('SECTORS:\n');
fprintf('  Informal/Formal wage ratio: %.4f\n\n', params.wI_wF_ratio);

fprintf('NUMERICAL:\n');
fprintf('  Asset grid points:         %d\n', params.N);
fprintf('  Asset range:               [%.2f, %.2f]\n\n', params.abar, params.amax);

%% Solve model
fprintf('Starting solution...\n\n');
tic;
results = aiyagari_2sectors_labor(params);
elapsed = toc;

fprintf('\nSolution time: %.2f seconds\n', elapsed);

%% Display results
fprintf('\n=================================================\n');
fprintf('              EQUILIBRIUM RESULTS\n');
fprintf('=================================================\n\n');

fprintf('PRICES:\n');
fprintf('  Interest rate:             %.4f (%.2f%%)\n', results.r, results.r*100);
fprintf('  Formal wage:               %.4f\n', results.w_formal);
fprintf('  Informal wage:             %.4f\n', results.w_informal);
fprintf('  Tax rate:                  %.4f (%.2f%%)\n\n', results.tau, results.tau*100);

fprintf('AGGREGATES:\n');
fprintf('  Capital:                   %.4f\n', results.K);
fprintf('  Formal labor (effective):  %.4f\n', results.L_formal);
fprintf('  Informal labor (effective):%.4f\n', results.L_informal);
fprintf('  Total labor (effective):   %.4f\n', results.L_total);
fprintf('  Informality rate (labor):  %.2f%%\n\n', ...
        100 * results.L_informal / results.L_total);

% Output
Y = results.K^params.alpha * results.L_total^(1-params.alpha);
fprintf('  Output:                    %.4f\n', Y);
fprintf('  Capital/Output ratio:      %.4f\n', results.K / Y);

% Consumption
C_total = results.g' * results.c;
fprintf('  Aggregate consumption:     %.4f\n', C_total);
fprintf('  Consumption/Output:        %.4f\n\n', C_total / Y);

%% Visualize results
figure('Position', [100 100 1200 800]);

% Panel 1: Value function
subplot(2,3,1);
plot(results.a, results.v(1:params.N), 'b-', 'LineWidth', 2); hold on;
plot(results.a, results.v(params.N+1:2*params.N), 'r--', 'LineWidth', 2);
xlabel('Assets (a)');
ylabel('Value function');
title('Value Function');
legend('Low prod (z_1)', 'High prod (z_2)', 'Location', 'best');
grid on;

% Panel 2: Consumption policy
subplot(2,3,2);
plot(results.a, results.c_z1, 'b-', 'LineWidth', 2); hold on;
plot(results.a, results.c_z2, 'r--', 'LineWidth', 2);
xlabel('Assets (a)');
ylabel('Consumption');
title('Consumption Policy');
legend('Low prod', 'High prod', 'Location', 'best');
grid on;

% Panel 3: Formal labor supply
subplot(2,3,3);
plot(results.a, results.lF_z1, 'b-', 'LineWidth', 2); hold on;
plot(results.a, results.lF_z2, 'r--', 'LineWidth', 2);
xlabel('Assets (a)');
ylabel('Formal labor (l_F)');
title('Formal Labor Supply');
legend('Low prod', 'High prod', 'Location', 'best');
grid on;

% Panel 4: Informal labor supply
subplot(2,3,4);
plot(results.a, results.lI_z1, 'b-', 'LineWidth', 2); hold on;
plot(results.a, results.lI_z2, 'r--', 'LineWidth', 2);
xlabel('Assets (a)');
ylabel('Informal labor (l_I)');
title('Informal Labor Supply');
legend('Low prod', 'High prod', 'Location', 'best');
grid on;

% Panel 5: Savings
subplot(2,3,5);
plot(results.a, results.s(1:params.N), 'b-', 'LineWidth', 2); hold on;
plot(results.a, results.s(params.N+1:2*params.N), 'r--', 'LineWidth', 2);
yline(0, 'k--');
xlabel('Assets (a)');
ylabel('Savings (s)');
title('Savings Policy');
legend('Low prod', 'High prod', 'Location', 'best');
grid on;

% Panel 6: Wealth distribution
subplot(2,3,6);
plot(results.a, results.g_z1, 'b-', 'LineWidth', 2); hold on;
plot(results.a, results.g_z2, 'r--', 'LineWidth', 2);
plot(results.a, results.g_z1 + results.g_z2, 'k-', 'LineWidth', 2);
xlabel('Assets (a)');
ylabel('Density');
title('Wealth Distribution');
legend('Low prod', 'High prod', 'Total', 'Location', 'best');
grid on;

sgtitle('Aiyagari Model with Two Sectors - Equilibrium');

%% Additional statistics
fprintf('=================================================\n');
fprintf('         DISTRIBUTION STATISTICS\n');
fprintf('=================================================\n\n');

% Wealth distribution
mean_wealth = results.g' * [results.a; results.a];
var_wealth = results.g' * ([results.a; results.a].^2) - mean_wealth^2;
std_wealth = sqrt(var_wealth);

fprintf('WEALTH:\n');
fprintf('  Mean:                      %.4f\n', mean_wealth);
fprintf('  Std deviation:             %.4f\n', std_wealth);
fprintf('  Coefficient of variation:  %.4f\n\n', std_wealth / mean_wealth);

% Labor supply distribution
mean_lF = results.g' * results.lF;
mean_lI = results.g' * results.lI;
mean_l_total = mean_lF + mean_lI;

fprintf('LABOR SUPPLY:\n');
fprintf('  Mean formal:               %.4f\n', mean_lF);
fprintf('  Mean informal:             %.4f\n', mean_lI);
fprintf('  Mean total:                %.4f\n', mean_l_total);
fprintf('  Informality (hours):       %.2f%%\n\n', 100 * mean_lI / mean_l_total);

% Income distribution
income_formal = results.w_formal * (1 - results.tau) * ...
                [results.z(1) * results.lF_z1; results.z(2) * results.lF_z2];
income_informal = results.w_informal * ...
                  [results.z(1) * results.lI_z1; results.z(2) * results.lI_z2];
income_total = income_formal + income_informal + results.r * [results.a; results.a];

mean_income = results.g' * income_total;
var_income = results.g' * (income_total.^2) - mean_income^2;
std_income = sqrt(var_income);

fprintf('INCOME:\n');
fprintf('  Mean:                      %.4f\n', mean_income);
fprintf('  Std deviation:             %.4f\n', std_income);
fprintf('  Coefficient of variation:  %.4f\n\n', std_income / mean_income);

%% Save results
save('aiyagari_2sectors_results.mat', 'results', 'params');
fprintf('\nResults saved to: aiyagari_2sectors_results.mat\n');

fprintf('\n=================================================\n');
fprintf('                   COMPLETE\n');
fprintf('=================================================\n');
