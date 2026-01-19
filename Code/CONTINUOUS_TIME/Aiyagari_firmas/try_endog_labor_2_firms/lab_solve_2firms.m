function eq = lab_solve_2firms(c, params)
% LAB_SOLVE_2FIRMS - Residual equation for 2-sector model (for use with fzero)
%
% Following Moll's lab_solve.m pattern exactly.
%
% The FOCs with SEPARABLE disutility are:
%   ell_F = (c^(-gamma) * w_F * z / psi_F)^phi
%   ell_I = (c^(-gamma) * w_I * theta * z / psi_I)^phi
%
% Budget constraint with zero savings:
%   c = w_F*z*ell_F + w_I*theta*z*ell_I + r*a
%
% Substituting FOCs into budget:
%   c = (w_F*z)^(1+phi) * c^(-gamma*phi) / psi_F^phi
%     + (w_I*theta*z)^(1+phi) * c^(-gamma*phi) / psi_I^phi
%     + r*a
%
% Residual: eq = 0 when c is the solution
%
% Inputs:
%   c      - consumption (scalar, the unknown)
%   params - [a, z, w_F, w_I, theta, r, gamma, phi, psi_F, psi_I]
%
% Output:
%   eq     - residual (should be zero at solution)

% Unpack parameters
a = params(1);
z = params(2);
w_F = params(3);
w_I = params(4);
theta = params(5);
r = params(6);
gamma = params(7);
phi = params(8);
psi_F = params(9);
psi_I = params(10);

% Coefficients from FOC substitution
A_F = (w_F * z)^(1 + phi) / psi_F^phi;
A_I = (w_I * theta * z)^(1 + phi) / psi_I^phi;

% Asset income
y_a = r * a;

% Residual equation: c - (A_F + A_I) * c^(-gamma*phi) - r*a = 0
eq = c - (A_F + A_I) * c^(-gamma * phi) - y_a;

end
