function eq = lab_solve_dual(l, params)
% LAB_SOLVE_DUAL - Solve dual problem for labor supply in two sectors
%
% INPUTS:
%   l      = [lf; li] = labor formal and informal (2x1 vector)
%   params = [a, z, wF, wI, r, gamma, frisch] (7x1 vector)
%            a      : assets
%            z      : productivity
%            wF     : formal wage
%            wI     : informal wage
%            r      : interest rate
%            gamma  : risk aversion
%            frisch : Frisch elasticity
%
% OUTPUTS:
%   eq = [eq1; eq2] = system of FOCs (2x1 vector)
%
% First Order Conditions:
%   eq1: lf^(1/φ) * c^γ = wF*z
%   eq2: li^(1/φ) * c^γ = wI*z
%
% where c = wF*z*lf + wI*z*li + r*a

% Extract parameters
a     = params(1);
z     = params(2);
wF    = params(3);
wI    = params(4);
r     = params(5);
gamma = params(6);
phi   = params(7);  % Frisch elasticity

% Extract labor supplies
lf = max(l(1), 0);  % Formal labor (non-negative)
li = max(l(2), 0);  % Informal labor (non-negative)

% Total income (consumption)
c = wF*z*lf + wI*z*li + r*a;

% Ensure consumption is positive
c = max(c, 1e-10);

% First order conditions
eq1 = (lf)^(1/phi) * c^gamma - wF*z;
eq2 = (li)^(1/phi) * c^gamma - wI*z;

eq = [eq1; eq2];

end