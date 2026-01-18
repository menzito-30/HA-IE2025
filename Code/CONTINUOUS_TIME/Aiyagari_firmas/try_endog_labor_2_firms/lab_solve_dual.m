function eq = lab_solve_dual(l, params)
% LAB_SOLVE_DUAL - Solve dual problem for labor supply in two sectors with taxes and transfers
%
% INPUTS:
%   l      = [lf; li] = labor formal and informal (2x1 vector)
%   params = [a, z, wF, wI, r, gamma, frisch, tau, tr] (9x1 vector)
%            a      : assets
%            z      : productivity
%            wF     : formal wage (pre-tax)
%            wI     : informal wage
%            r      : interest rate
%            gamma  : risk aversion
%            frisch : Frisch elasticity
%            tau    : tax rate on formal
%            tr     : lump-sum transfer
%
% OUTPUTS:
%   eq = [eq1; eq2] = system of FOCs (2x1 vector)
%
% First Order Conditions (intensive margin):
%   eq1: lf^(1/ϕ) * c^γ = (1 - τ) wF z
%   eq2: li^(1/ϕ) * c^γ = wI z
%
% where c = (1 - τ) wF z lf + wI z li + r a + tr

% Extract parameters
a     = params(1);
z     = params(2);
wF    = params(3);
wI    = params(4);
r     = params(5);
gamma = params(6);
phi   = params(7);  % Frisch elasticity
tau   = params(8);
tr    = params(9);

% Extract labor supplies (non-negative)
lf = max(l(1), 0);  % Formal labor
li = max(l(2), 0);  % Informal labor

% Total income (consumption, after-tax + transfer)
c = (1 - tau) * wF * z * lf + wI * z * li + r * a + tr;

% Ensure consumption is positive
c = max(c, 1e-10);

% First order conditions (adjusted for tax)
eq1 = (lf) ^ (1 / phi) * c ^ gamma - (1 - tau) * wF * z;
eq2 = (li) ^ (1 / phi) * c ^ gamma - wI * z;

eq = [eq1; eq2];

end