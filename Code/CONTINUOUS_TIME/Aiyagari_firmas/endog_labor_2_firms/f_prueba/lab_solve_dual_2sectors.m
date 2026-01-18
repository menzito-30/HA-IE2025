function [lF, lI, c, converged] = lab_solve_dual_2sectors(a, z, wF, wI, tau, r, gamma, frisch, vpa)
% LAB_SOLVE_DUAL_2SECTORS - Solve dual problem for labor supply in formal and informal sectors
%
% INPUTS:
%   a      : assets
%   z      : productivity
%   wF     : formal wage
%   wI     : informal wage  
%   tau    : proportional tax on formal income
%   r      : interest rate
%   gamma  : risk aversion
%   frisch : Frisch elasticity
%   vpa    : marginal value of assets v'(a)
%
% OUTPUTS:
%   lF        : formal labor supply
%   lI        : informal labor supply
%   c         : consumption
%   converged : convergence flag
%
% FOCs (intensive margin):
%   c^(-gamma) = v'(a)
%   lF^(1/frisch) * c^gamma = wF*(1-tau)*z
%   lI^(1/frisch) * c^gamma = wI*z

% Initialize
converged = true;
phi = frisch;

% Case 1: Try interior solution (both lF > 0 and lI > 0)
wF_net = wF * (1 - tau) * z;  % After-tax formal wage
wI_net = wI * z;              % Informal wage (no tax)

% Initial guess
l0 = [0.3; 0.3];  % [lF; lI]

% Solve system of FOCs
options = optimset('Display', 'off', 'TolX', 1e-10, 'TolFun', 1e-10);
try
    [l_sol, ~, exitflag] = fsolve(@(l) foc_system(l, a, wF_net, wI_net, r, gamma, phi), ...
                                   l0, options);
    
    if exitflag > 0 && all(l_sol > 1e-6)
        % Interior solution found
        lF = l_sol(1);
        lI = l_sol(2);
        c = wF_net * lF + wI_net * lI + r * a;
        return;
    end
catch
    % Fall through to corner solutions
end

% Case 2: Only formal (lF > 0, lI = 0)
try
    lF = fzero(@(l) foc_formal_only(l, a, wF_net, r, gamma, phi), 0.5, options);
    if lF > 1e-6
        lI = 0;
        c = wF_net * lF + r * a;
        if c > 1e-6
            return;
        end
    end
catch
    % Fall through
end

% Case 3: Only informal (lF = 0, lI > 0)
try
    lI = fzero(@(l) foc_informal_only(l, a, wI_net, r, gamma, phi), 0.5, options);
    if lI > 1e-6
        lF = 0;
        c = wI_net * lI + r * a;
        if c > 1e-6
            return;
        end
    end
catch
    % Fall through
end

% Case 4: No work (corner at borrowing constraint)
lF = 0;
lI = 0;
c = max(r * a, 1e-10);
converged = false;

end

%% Nested functions for FOC systems

function eq = foc_system(l, a, wF_net, wI_net, r, gamma, phi)
    % System for interior solution
    lF = max(l(1), 0);
    lI = max(l(2), 0);
    
    c = wF_net * lF + wI_net * lI + r * a;
    c = max(c, 1e-10);
    
    % FOCs
    eq1 = (lF)^(1/phi) * c^gamma - wF_net;
    eq2 = (lI)^(1/phi) * c^gamma - wI_net;
    
    eq = [eq1; eq2];
end

function eq = foc_formal_only(lF, a, wF_net, r, gamma, phi)
    % FOC when only working in formal sector
    lF = max(lF, 0);
    c = wF_net * lF + r * a;
    c = max(c, 1e-10);
    
    eq = (lF)^(1/phi) * c^gamma - wF_net;
end

function eq = foc_informal_only(lI, a, wI_net, r, gamma, phi)
    % FOC when only working in informal sector
    lI = max(lI, 0);
    c = wI_net * lI + r * a;
    c = max(c, 1e-10);
    
    eq = (lI)^(1/phi) * c^gamma - wI_net;
end
