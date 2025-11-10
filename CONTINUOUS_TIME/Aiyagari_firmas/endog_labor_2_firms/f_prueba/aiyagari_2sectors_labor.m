function results = aiyagari_2sectors_labor(params)
% AIYAGARI_2SECTORS_LABOR - Solve Aiyagari model with endogenous labor supply
%                           in formal and informal sectors
%
% Based on Pontus Rendahl's continuous time methods and Benjamin Moll's
% labor supply implementation
%
% STRUCTURE:
%   - 2 sectors: Formal (taxed) and Informal (untaxed)
%   - Endogenous labor supply (intensive margin)
%   - Idiosyncratic productivity shocks (2-state Poisson)
%   - Government with balanced budget
%
% INPUTS:
%   params : structure with model parameters
%
% OUTPUTS:
%   results : structure with equilibrium objects

%% Extract parameters
rho     = params.rho;        % Discount rate
gamma   = params.gamma;      % Risk aversion
frisch  = params.frisch;     % Frisch elasticity
alpha   = params.alpha;      % Capital share
delta   = params.delta;      % Depreciation
mu      = params.mu;         % UI replacement rate (not used here)
lambda1 = params.lambda1;    % Transition intensity z1 -> z2
lambda2 = params.lambda2;    % Transition intensity z2 -> z1
z1      = params.z1;         % Low productivity
z2      = params.z2;         % High productivity
abar    = params.abar;       % Borrowing constraint
N       = params.N;          % Grid points for assets
amax    = params.amax;       % Max assets
Gamma   = params.Gamma;      % Implicit method parameter

wI_wF_ratio = params.wI_wF_ratio;  % Wage ratio informal/formal
tau_target  = params.tau_target;   % Target tax rate (or empty for endogenous)

%% Setup grids
a = linspace(abar, amax, N)';
da = a(2) - a(1);

z = [z1; z2];  % Productivity states
Nz = length(z);

%% Build difference operators (sparse for efficiency)
% Forward difference
Df = spdiags([-ones(N,1), ones(N,1)], 0:1, N, N) / da;
Df(N,N-1:N) = [0, -1/da];  % Adjust last row

% Backward difference  
Db = spdiags([-ones(N,1), ones(N,1)], -1:0, N, N) / da;
Db(1,1:2) = [1/da, 0];     % Adjust first row

%% Build productivity transition matrix B (2N × 2N)
B = sparse(2*N, 2*N);
B(1:N, 1:N) = -lambda1 * speye(N);
B(1:N, N+1:2*N) = lambda1 * speye(N);
B(N+1:2*N, 1:N) = lambda2 * speye(N);
B(N+1:2*N, N+1:2*N) = -lambda2 * speye(N);

%% Find stationary productivity distribution
s_prod = zeros(2,1);
s_prod(1) = lambda2 / (lambda1 + lambda2);  % Prob of z1
s_prod(2) = lambda1 / (lambda1 + lambda2);  % Prob of z2

fprintf('\n=== Solving Aiyagari Model with Two Sectors ===\n');
fprintf('Stationary productivity distribution: [%.4f, %.4f]\n', s_prod(1), s_prod(2));

%% Outer loop: Find equilibrium interest rate
r_min = -delta + 0.001;
r_max = rho - 0.001;
r_guess = 0.03;

max_iter_r = 50;
tol_r = 1e-4;

for iter_r = 1:max_iter_r
    
    r = r_guess;
    
    % Implied wage from firm FOC
    % w = (1-alpha) * (r+delta/alpha)^(alpha/(alpha-1))
    w_formal = (1 - alpha) * ((r + delta) / alpha)^(alpha / (alpha - 1));
    w_informal = wI_wF_ratio * w_formal;
    
    % Tax rate (exogenous for now)
    if isempty(tau_target)
        tau = 0.15;  % Default
    else
        tau = tau_target;
    end
    
    fprintf('\n--- Iteration %d: r = %.4f, wF = %.4f, wI = %.4f, tau = %.4f ---\n', ...
            iter_r, r, w_formal, w_informal, tau);
    
    %% Solve HJB equation using implicit method
    
    % Initial guess for value function
    v0 = zeros(2*N, 1);
    c_guess = r * a + 0.1;
    for iz = 1:Nz
        v0((iz-1)*N+1:iz*N) = (c_guess.^(1-gamma)) / (1-gamma);
    end
    
    % Iteration parameters
    max_iter_v = 1000;
    tol_v = 1e-6;
    
    for iter_v = 1:max_iter_v
        
        v = v0;
        
        % Storage for policy functions and drift
        c = zeros(2*N, 1);
        lF = zeros(2*N, 1);
        lI = zeros(2*N, 1);
        s = zeros(2*N, 1);
        
        % Loop over productivity states
        for iz = 1:Nz
            
            idx = (iz-1)*N+1:iz*N;
            v_z = v(idx);
            
            % Compute v'(a) using forward and backward differences
            vpa_f = Df * v_z;
            vpa_b = Db * v_z;
            
            % Boundary conditions
            vpa_b(1) = (w_formal * (1-tau) * z(iz) * 0.5 + r * abar)^(-gamma);
            vpa_f(N) = (w_formal * (1-tau) * z(iz) * 0.5 + r * amax)^(-gamma);
            
            % Loop over asset grid
            for ia = 1:N
                
                % Solve for optimal labor and consumption
                [lF_f, lI_f, c_f, ~] = lab_solve_dual_2sectors(a(ia), z(iz), ...
                    w_formal, w_informal, tau, r, gamma, frisch, vpa_f(ia));
                
                [lF_b, lI_b, c_b, ~] = lab_solve_dual_2sectors(a(ia), z(iz), ...
                    w_formal, w_informal, tau, r, gamma, frisch, vpa_b(ia));
                
                % Compute savings (drift)
                s_f = w_formal * (1-tau) * z(iz) * lF_f + w_informal * z(iz) * lI_f + r * a(ia) - c_f;
                s_b = w_formal * (1-tau) * z(iz) * lF_b + w_informal * z(iz) * lI_b + r * a(ia) - c_b;
                
                % Upwind scheme
                if s_f > 0
                    c(idx(ia)) = c_f;
                    lF(idx(ia)) = lF_f;
                    lI(idx(ia)) = lI_f;
                    s(idx(ia)) = s_f;
                elseif s_b < 0
                    c(idx(ia)) = c_b;
                    lF(idx(ia)) = lF_b;
                    lI(idx(ia)) = lI_b;
                    s(idx(ia)) = s_b;
                else
                    % No drift at this point
                    c(idx(ia)) = r * a(ia) + w_formal * (1-tau) * z(iz) * 0.3 + ...
                                 w_informal * z(iz) * 0.3;
                    lF(idx(ia)) = 0.3;
                    lI(idx(ia)) = 0.3;
                    s(idx(ia)) = 0;
                end
                
            end
        end
        
        % Build matrix SD
        SD = spdiags(s, 0, 2*N, 2*N);
        
        % Build drift matrix for each productivity state
        P = sparse(2*N, 2*N);
        for iz = 1:Nz
            idx = (iz-1)*N+1:iz*N;
            s_z = s(idx);
            
            % Indicators
            If = (s_z > 0);
            Ib = (s_z < 0);
            
            % Upwind difference operator
            D_upwind = spdiags(If, 0, N, N) * Df + spdiags(Ib, 0, N, N) * Db;
            
            % Drift matrix for this state
            SD_z = spdiags(s_z, 0, N, N) * D_upwind;
            
            P(idx, idx) = SD_z;
        end
        
        % Add productivity transitions
        P = P + B;
        
        % Utility flow
        l_total = lF + lI;
        u_flow = (c.^(1-gamma)) / (1-gamma) - (l_total.^(1 + 1/frisch)) / (1 + 1/frisch);
        
        % Update value function (implicit method)
        A = (rho + 1/Gamma) * speye(2*N) - P;
        b = u_flow + v / Gamma;
        
        v_new = A \ b;
        
        % Check convergence
        error_v = max(abs(v_new - v));
        
        if mod(iter_v, 100) == 0 || iter_v == 1
            fprintf('  HJB iter %4d: error = %.2e\n', iter_v, error_v);
        end
        
        if error_v < tol_v
            fprintf('  HJB converged in %d iterations\n', iter_v);
            break;
        end
        
        v0 = v_new;
        
    end
    
    if iter_v == max_iter_v
        warning('HJB did not converge');
    end
    
    %% Solve for stationary distribution (Kolmogorov Forward Equation)
    
    % 0 = P' * g
    % Use trick: replace first row with normalization
    PT = P';
    PT_hat = PT;
    PT_hat(1,:) = ones(1, 2*N);  % Normalization: sum(g) = 1
    
    b_g = zeros(2*N, 1);
    b_g(1) = 1;
    
    g = PT_hat \ b_g;
    g = g / sum(g);  % Ensure normalization
    
    % Check non-negativity
    if any(g < 0)
        warning('Negative densities detected');
        g = max(g, 0);
        g = g / sum(g);
    end
    
    %% Compute aggregate capital supply
    g_z1 = g(1:N);
    g_z2 = g(N+1:2*N);
    
    K_supply = g_z1' * a + g_z2' * a;
    
    % Aggregate labor (effective)
    L_formal = g_z1' * (z1 * lF(1:N)) + g_z2' * (z2 * lF(N+1:2*N));
    L_informal = g_z1' * (z1 * lI(1:N)) + g_z2' * (z2 * lI(N+1:2*N));
    L_total = L_formal + L_informal;
    
    fprintf('  Aggregates: K = %.4f, L_F = %.4f, L_I = %.4f, L_total = %.4f\n', ...
            K_supply, L_formal, L_informal, L_total);
    
    %% Compute implied capital demand from firm FOC
    % r = alpha * (K/L)^(alpha-1) - delta
    % => K/L = ((r + delta)/alpha)^(1/(alpha-1))
    
    K_over_L = ((r + delta) / alpha)^(1 / (alpha - 1));
    K_demand = K_over_L * L_total;
    
    fprintf('  Capital: Supply = %.4f, Demand = %.4f, Excess = %.4f\n', ...
            K_supply, K_demand, K_supply - K_demand);
    
    %% Check market clearing
    excess_K = K_supply - K_demand;
    
    if abs(excess_K) < tol_r
        fprintf('\n=== EQUILIBRIUM FOUND ===\n');
        fprintf('Interest rate: r = %.4f\n', r);
        fprintf('Capital: K = %.4f\n', K_supply);
        fprintf('Formal wage: wF = %.4f\n', w_formal);
        fprintf('Informal wage: wI = %.4f\n', w_informal);
        fprintf('Tax rate: tau = %.4f\n', tau);
        break;
    end
    
    %% Update interest rate (bisection)
    if excess_K > 0
        % Too much saving, lower r
        r_max = r;
    else
        % Too little saving, raise r
        r_min = r;
    end
    
    r_guess = (r_min + r_max) / 2;
    
end

if iter_r == max_iter_r
    warning('Outer loop did not converge');
end

%% Package results
results.r = r;
results.w_formal = w_formal;
results.w_informal = w_informal;
results.tau = tau;
results.K = K_supply;
results.L_formal = L_formal;
results.L_informal = L_informal;
results.L_total = L_total;

results.a = a;
results.z = z;
results.v = v;
results.c = c;
results.lF = lF;
results.lI = lI;
results.s = s;
results.g = g;

results.c_z1 = c(1:N);
results.c_z2 = c(N+1:2*N);
results.lF_z1 = lF(1:N);
results.lF_z2 = lF(N+1:2*N);
results.lI_z1 = lI(1:N);
results.lI_z2 = lI(N+1:2*N);
results.g_z1 = g_z1;
results.g_z2 = g_z2;

fprintf('\n=== Solution Complete ===\n');

end
