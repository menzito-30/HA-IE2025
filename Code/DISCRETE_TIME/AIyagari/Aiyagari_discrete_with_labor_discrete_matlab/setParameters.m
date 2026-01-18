%% Economic Parameters
% ---------------------------------------
% Household
par.beta        = 0.99;     % Discount factor - calibrated to match SS interest rate
par.sigma       = 1;        % CRRA
par.psi1        = 1;        % disutility of work
par.psi2        = 2;        % Inverse Frisch elasticity

% Income Process
par.rhoz        = 0.96566;        % Persistence of productivity
par.sigmaz      = sqrt(0.01695);  % STD of productivity shocks
par.uncond_muz  = 0;              % unconditional mean of z for rouwenhorst
par.uncond_sigmaz   = par.sigmaz/sqrt(1-par.rhoz^2);    % unconditional stddev of z for rouwenhorst

% Firms
par.alpha       = 0.36;     % capital share
par.delta       = 0.025;    % capital depreciation

%% Grid Parameters
% ---------------------------------------
% Idiosyncratic States
mpar.na         = 100;      % grid points individual capital
mpar.nz         = 7;        % grid points individual productivity

%% Numerical Parameters
% ---------------------------------------
mpar.crit    = 1e-10;

