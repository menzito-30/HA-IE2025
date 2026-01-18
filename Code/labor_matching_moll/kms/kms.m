%==========================================================================
%             Statonary Equilibrium in Krusell-Mukoyama-Sahin
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

% Grid for assets
gridd.I    = 1001;
gridd.amin = 0;
gridd.amax = 150;

% Numerical parameters
num.tol_hjb   = 1e-8;  % tolerance level for HJB equations
num.tol_ent   = 1e-6;  % tolerance level for free entry condition
num.tol_mkt   = 1e-6;  % tolerance level for asset market clearing
num.maxit_hjb = 50;    % maximum number of iterations when solving HJB equation
num.maxit_out = 5000;  % maximum number of iterations in outer loop
num.Delta     = 2000;  % step size in HJB equation
num.dtheta    = 1e-1;  % step size in updating theta
num.dk        = 1e-3;  % step size in updating k

% Initial guesses
guess.theta = 1;
guess.k     = 25.61;

% Bargaining solution
% bargain  = 'nash';
% param.xi = 0.395; % with Nash

bargain  = 'egalitarian';
param.xi = 0.199; % with egalitarian

%=========================================================================
%%                         Solve for Equilibrium
%=========================================================================
tic
st = kmsfun(param, gridd, num, guess, bargain);
toc

%=========================================================================
%%                              Graphs
%=========================================================================
a    = linspace(gridd.amin, gridd.amax, gridd.I);
acut = floor(gridd.I/1.5);

% WORKER VALUE FUNCTION
figure(1)
plot(a, st.W, 'linewidth', 1.5)
box on; grid on; xlim([-5, a(acut)])
title('Worker value functions'); xlabel('assets'); ylabel('value function')
legend('employed', 'unemployed', 'Location', 'northwest');

% FILLED JOB VALUE FUNCTION
figure(2)
plot(a, st.J, 'linewidth', 1.5)
box on; grid on; xlim([-5, a(acut)])
title('Filled job value function'); xlabel('assets'); ylabel('value function')

% CONSUMPTION
figure(3)
hold on; plot(a, st.c, 'linewidth', 1.5); hold off;
box on; grid on; xlim([-5, a(acut)])
title('Consumption'); xlabel('assets'); ylabel('consumption')
legend('employed', 'unemployed', 'location', 'southeast');

% SAVINGS - can be used to check whether amax is high enough (negative for employed)
figure(4)
hold on; plot(a, st.adot, 'linewidth', 1.5); plot(a, zeros(length(a), 1), 'k--'); hold off;
lgnd = legend('employed', 'unemployed', 'location', 'best');
box on; grid on; xlim([-5, a(acut)])
title('Savings'); xlabel('assets'); ylabel('savings')

% STATIONARY DISTRIBUTION (marginal)
figure(5)
gg_marg       = zeros(gridd.I, 2);
gg_marg(:, 1) = st.gg(:, 1) / (1-st.u);
gg_marg(:, 2) = st.gg(:, 2) / st.u;
hold on; plot(a, gg_marg, 'linewidth', 1.5); hold off;
legend('employed', 'unemployed','location', 'northeast');
box on; grid on; xlim([-5, a(acut)])
title('Marginals of the stationary distribution'); xlabel('assets'); ylabel('density')

% WAGE
figure(6)
plot(a, st.w, 'linewidth', 1.5);
box on; grid on; xlim([-5, a(acut)]);
title('Wage schedule'); xlabel('assets'); ylabel('wage');