function [pol, JD, SS] = main_steadystate(par, mpar, grid, mesh, Pz)
% _________________________________________________________________________
%% Solve for steady state
% _________________________________________________________________________

% 0) solve for aggregate per-capita-capital
% ---------------------------------------
tic
f = @(KN) compute_steady_state(KN, par, mpar, grid, mesh, Pz);
options = optimoptions(@fsolve,'Display','iter');
KN = fsolve(f,grid.KN,options);
grid.KN = KN;
toc

% 1) Compute factor returns
% ---------------------------------------
[SS.R, SS.W] = factor_returns(KN, par, grid);

% 2) Guess initial policy function
% ---------------------------------------
c_guess = policyguess(mpar, mesh, SS.R, SS.W);

% 3) Solve for policy functions
% ---------------------------------------
pol = policies_EGM(par, mpar, grid, mesh, c_guess, Pz, SS.R, SS.W);

% 4) Solve for stationary distribution
% ---------------------------------------
JD = JDiteration(mpar, grid, Pz, pol);

end



% _________________________________________________________________________
%% Function: compute_steady_state solves for SS policy functions and distribution
% _________________________________________________________________________
function resid = compute_steady_state(KN, par, mpar, grid, mesh, Pz)

% 1) Compute factor returns (prices in HH BC)
% ---------------------------------------
[R, W] = factor_returns(KN, par, grid);

% 2) Guess initial policy function
% ---------------------------------------
c_guess = policyguess(mpar, mesh, R, W);

% 3) Solve for policy functions
% ---------------------------------------
pol = policies_EGM(par, mpar, grid, mesh, c_guess, Pz, R, W);

% 4) Solve for stationary distribution
% ---------------------------------------
JD = JDiteration(mpar, grid, Pz, pol);

% 5) Compute residual capital-to-labor
% ---------------------------------------
AggCapital = grid.a*sum(JD,2);
AggLabor = sum(sum(mesh.z .* pol.l .* JD));
resid = AggCapital/AggLabor - grid.KN;

end