%==========================================================================
%            Transition Dynamics in Krusell-Mukoyama-Sahin
%
%
% Author: Bence Bardoczy
% Date  : 10/12/2017
%==========================================================================

%=========================================================================
%%                           Housekeeping
%=========================================================================
close all;
clear all;
clc;

%=========================================================================
%%                             Initialization
%=========================================================================
% Structural parameters, unit of time is quarter
param.gam   = 1;      % relative risk aversion
param.rho   = 0.01;   % discount rate
param.chi   = 1.7935; % matching efficiency
param.sig   = 0.1038; % separation rate
param.eta   = 0.72;   % matching elasticity to unemp.
param.bet   = 0.72;   % worker bargaining power
param.h     = 0.75;   % home production
param.alfa  = 0.3;    % capital share
param.delta = 0.021;  % depreciation rate
param.z     = 1;      % steady state TFP

% Grid: asset
gridd.I     = 501;
gridd.amin  = 0;
gridd.amax  = 150;

% Bargaining solution
% bargain  = 'nash';
% param.xi = 0.395; % with Nash

bargain  = 'egalitarian';
param.xi = 0.199; % with egalitarian

%=========================================================================
%%                   Solve for Stationary Equilibrium
%=========================================================================
% Numerical parameters for steady state calculation
num1.tol_hjb   = 1e-8; % tolerance level for HJB equations
num1.tol_ent   = 1e-6;  % tolerance level for free entry condition
num1.tol_mkt   = 1e-6;  % tolerance level for asset market clearing
num1.maxit_hjb = 50;    % maximum number of iterations when solving HJB equation
num1.maxit_out = 5000;  % maximum number of iterations in outer loop
num1.Delta     = 2000;  % step size in HJB equation
num1.dtheta    = 1e-1;  % step size in updating theta
num1.dk        = 1e-3;  % step size in updating r

% Initial guesses for steady state calculation
guess1.theta = 1;
guess1.k     = 1.2 * (param.alfa / (param.rho+param.delta)) ^ (1/(1-param.alfa));

% Call function to do the iteration
tic
st = kmsfun(param, gridd, num1, guess1, bargain);
toc

%=========================================================================
%%                   Solve for Transition Dynamics
%=========================================================================
% TFP Shock
shock.cor = 0.5;
shock.z0  = 1.02;

% Grid: time,
% time steps are small in the beginning to capture the action right after the shock, then get longer to gain speed
gridd.d0    = 1/3;
gridd.d1    = 5;
gridd.dtvec = [ones(60, 1)*gridd.d0; linspace(gridd.d0, gridd.d1, 60)']; % time steps
gridd.time  = [0; cumsum(gridd.dtvec)];       % non-uniform time grid
gridd.T     = gridd.time(end);                % horizon in quarters
gridd.N     = length(gridd.time);             % number of grid points
% plot(1:gridd.N, gridd.time); xlabel('number of grid points'); ylabel('quarters covered')

% Numerical parameters for steady state calculation
num2.tol_ent   = 1e-5;  % tolerance level for free entry condition
num2.tol_mkt   = 1e-5;  % tolerance level for asset market clearing
num2.tol_eqty  = 1e-5;  % tolerance level for equity pric on impact
num2.maxit_out = 10000; % maximum number of iterations in outer loop
num2.dtheta    = 1e-1 * exp(-linspace(0, 3, gridd.N))';        % step size for updating theta_t
num2.dk        = [0; 1e-2 * exp(-linspace(0, 3, gridd.N-1))']; % step size for updating k_t
num2.dp0       = 1;                                            % step size for updating p0

% Initial guesses for steady state calculation
guess2.theta_t = ones(gridd.N, 1) * st.theta;
guess2.k_t     = ones(gridd.N, 1) * st.k;
guess2.p0      = 0.9 * st.p;
% load guess_k
% guess2.theta_t = theta_t;
% guess2.k_t     = k_t;
% guess2.p0      = p0;
% clearvars theta_t k_t p0


% Call function to do the iteration
tic
trans = kmstranfun(param, gridd, num2, guess2, shock, st, bargain);
toc

%=========================================================================
%%                          Impulse responses
%=========================================================================
N    = gridd.N;
da   = (gridd.amax-gridd.amin) / (gridd.I-1);
a    = linspace(gridd.amin, gridd.amax, gridd.I);
acut = floor(gridd.I*0.5);
tcut = floor(N*0.1);
time = [linspace(-5, -0.5, 21)'; gridd.time];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LABOR MARKET
figure(1)

% unemployment
subplot(2, 2, 2)
irf.u =  [ones(21, 1)*st.u; trans.u_t];
plot(time, irf.u, 'linewidth', 1.5); hold on; plot(time, st.u*ones(length(time), 1), 'r--'); hold off;
box on; grid on; xlim([-5, 20]); %ylim([-6, 1])
title('Unemployment rate'); xlabel('quarters'); ylabel('%')

% job-finding rate
subplot(2, 2, 1)
st.fmonth      = 1 - exp(-st.f/3);
trans.fmonth_t = 1 - exp(-trans.f_t/3);
irf.f = [ones(21, 1)*st.fmonth; trans.fmonth_t];
plot(time, irf.f, 'linewidth', 1.5); hold on; plot(time, st.fmonth*ones(length(time), 1), 'r--'); hold off;
box on; grid on; xlim([-5, 20]); %ylim([0.445, 0.48])
title('Job-finding rate'); xlabel('quarters'); ylabel('monthly prob.')

% wage on impact
subplot(2, 2, 3)
hold on
plot(a, trans.w_t{1}, 'linewidth', 1.5)
plot(a, st.w, 'r--')
% plot(a, trans.w_t{60}, 'k-', 'linewidth', 1.5)
hold off
box on; grid on; xlim([-5, 70]); %ylim([1.78, 1.87])
% legend('stationary', 'on impact', 'Location', 'southeast')
title('Wage schedule on impact'); xlabel('assets'); ylabel('w')

% average wage
subplot(2, 2, 4)
trans.meanwage_t = zeros(N, 1);
for n = 1:N
   trans.meanwage_t(n) = (trans.gg_t{n}(1:gridd.I)' * trans.w_t{n}) / (1 - trans.u_t(n)) * da;
end
irf.meanwage =  100 * [zeros(21, 1); (trans.meanwage_t - st.meanwage)/st.meanwage];
plot(time, irf.meanwage, 'linewidth', 1.5); hold on; plot(time, zeros(length(time), 1), 'r--'); hold off;
box on; grid on; xlim([-5, 20]); %ylim([-0.1, 1])
title('Average wage bill'); xlabel('quarters'); ylabel('% deviation from ss')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ASSETS
figure(2)

% capital stock
subplot(2, 2, 1)
irf.K =  100 * [zeros(21, 1); (trans.K_t - st.K)/st.K];
plot(time, irf.K, 'linewidth', 1.5); hold on; plot(time, zeros(length(time), 1), 'r--'); hold off;
box on; grid on; xlim([-5, 20]); %ylim([-0.05, 0.38])
title('Capital'); xlabel('quarters'); ylabel('% deviation from ss')

% equity price
subplot(2, 2, 2)
irf.p =  100 * [zeros(21, 1); (trans.p_t - st.p)/st.p];
plot(time, irf.p, 'linewidth', 1.5); hold on; plot(time, zeros(length(time), 1), 'r--'); hold off;
box on; grid on; xlim([-5, 20]); %ylim([-5, 30])
title('Equity price'); xlabel('quarters'); ylabel('% deviation from ss')

% total assets
subplot(2, 2, 3)
A_t = trans.K_t + trans.p_t;
irf.A = 100 * [zeros(21, 1); (A_t - st.K - st.p)/(st.K+st.p)];
plot(time, irf.A, 'linewidth', 1.5); hold on; plot(time, zeros(length(time), 1), 'r--'); hold off;
box on; grid on; xlim([-5, 20]); %ylim([-0.05, 0.38])
title('Total assets'); xlabel('quarters'); ylabel('% deviation from ss')

% interest rate
subplot(2, 2, 4)
st.rannual       = 1 - exp(-4*(st.r-param.delta));
trans.rannual_t = 1 - exp(-4*(trans.r_t-param.delta));
irf.r =  100 * [st.rannual * ones(21, 1); trans.rannual_t];
plot(time, irf.r, 'linewidth', 1.5); hold on; plot(time, 100*st.rannual*ones(length(time), 1), 'r--'); hold off;
box on; grid on; xlim([-5, 20]); %ylim([3.85, 4.15])
title('Interest rate'); xlabel('quarters'); ylabel('%, annualized')