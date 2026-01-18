function sol = kmstranfun(param, gridd, num, guess, shock, steady, bargain)
% Function to calculate the stationary equilibrium of a Krusell-Mukoyama-Sahin model.
%
% Input:
%     param  : struct, contains structural parameters of the model
%     gridd  : struct, specifies grid for assets and time
%     num    : struct, contains numerical parameters
%     guess  : struct, initial guess for theta_t, k_t
%     shock  : struct, specifies TFP shock
%     steady : struct, contains stationary equilibrium (initial and terminal conditions)
%
% Output:
%     sol : struct, model objects along transition

%==========================================================================
%                              Initialization
%==========================================================================
% Open up structures (is this slow?)
gam       = param.gam;
rho       = param.rho;
chi       = param.chi;
sig       = param.sig;
eta       = param.eta;
bet       = param.bet;
h         = param.h;
alfa      = param.alfa;
delta     = param.delta;
xi        = param.xi;
z         = param.z;

I         = gridd.I;
amin      = gridd.amin;
amax      = gridd.amax;
N         = gridd.N;
dtvec     = gridd.dtvec;
time      = gridd.time;

cor       = shock.cor;
z0        = shock.z0;

tol_ent   = num.tol_ent;
tol_mkt   = num.tol_mkt;
tol_eqty  = num.tol_eqty;
maxit_out = num.maxit_out;
dtheta    = num.dtheta;
dk        = num.dk;
dp0       = num.dp0;

theta_t   = guess.theta_t;
k_t       = guess.k_t;
k_t(1)    = steady.k; % capital is a state, cannot jump on impact
p0        = guess.p0;

% Finish building the asset grid
a    = linspace(amin, amax, I)';
aa   = [a, a];
da   = (amax-amin) / (I-1);

% Construct TFP path
nu  = 1 - cor;
z_t = (z0 - 1) * exp(-nu*time) + 1;

% Preallocation
FE_t = zeros(N, 1);
AD_t = zeros(N, 1);
u_t  = zeros(N, 1);
K_t  = zeros(N, 1);
d_t  = zeros(N, 1);
p_t  = zeros(N, 1);
w_t  = cell(N, 1);
W_t  = cell(N, 1);
J_t  = cell(N, 1);
gg_t = cell(N, 1);
Lambda_t = cell(N, 1);
profit_t = cell(N, 1);

%=========================================================================
%%                             Iteration
%=========================================================================
% Solve for transition given current guess for r_t, theta_t
tic
for it = 1:maxit_out
   disp(['Outer loop, Iteration ' int2str(it)]);

   %               1. Prices and labor market variables
   %---------------------------------------------------------------------
   % transition probabilities
   f_t = chi * theta_t.^(1-eta);
   q_t = chi * theta_t.^(-eta);

   % transition matrix
   for n = 1:N
      lam         = [-sig, sig; f_t(n), -f_t(n)];
      Lambda_t{n} = [speye(I)*lam(1,1), speye(I)*lam(1,2);
                     speye(I)*lam(2,1), speye(I)*lam(2,2)];
   end

   % prices
   r_t = alfa * z_t .* k_t .^ (alfa-1);

   %      2. Iterate HJB equations of household and firm backwards
   %---------------------------------------------------------------------
   % Start from terminal distribution
   W = steady.W;
   J = steady.J;
   w = steady.w;
   for n = N:-1:1
      % current variables
      k = k_t(n);
      r = r_t(n);
      z = z_t(n);
      Lambda = Lambda_t{n};

      % current stepsize (we have one more period)
      if n == 1
         dt = dtvec(1);
      else
         dt = dtvec(n-1);
      end

      % update functions of wage
      y = [w, ones(I, 1)*h];
      profit = z * k^alfa - r*k - w;

      %% A) WORKER
      %----------------------------------------------
      % interior partials
      WB(2:I, :)   = (W(2:I, :) - W(1:I-1, :)) / da;
      WF(1:I-1, :) = (W(2:I, :) - W(1:I-1, :)) / da;

      % state constraint boundary condtitions
      WB(1, :) = (y(1, :) + (r-delta)*amin).^(-gam);
      WF(I, :) = (y(I, :) + (r-delta)*amax).^(-gam);

      % consumption
      cB = max(WB, 1e-6).^(-1/gam);
      cF = max(WF, 1e-6).^(-1/gam);
      c0 = y + (r-delta)*aa;

      % savings
      sB = y + (r-delta)*aa - cB;
      sF = y + (r-delta)*aa - cF;
      sF(I, :) = 0; % make sure never will be used

      % cases (uses concavity)
      Ib = sB < 0;
      If = sF > 0;
      I0 = 1 - Ib - If;

      % upwinding
      c = If.*cF + Ib.*cB + I0.*c0;
      if gam == 1
         util = log(c);
      else
         util = c.^(1-gam) / (1-gam);
      end
      adot = y + (r-delta)*aa - c;

      % diagonals
      X = -min(adot, 0) / da;
      Z =  max(adot, 0) / da;
      Y = -X - Z;

      % prepare for spdiags
      X = X(2:I,:);
      Z = [zeros(1,2); Z(1:I-1,:)];

      % transition matrix
      for j = 1:2
         Ak{j} = spdiags(Y(:,j), 0, I, I) + spdiags(X(:,j), -1, I, I) + spdiags(Z(:,j), 1, I, I);
      end
      A = [Ak{1}, sparse(I,I); sparse(I,I), Ak{2}];

      % set up and solve linear system
      B = (rho + 1/dt)*speye(I*2) - A - Lambda;
      b = vec(util) + vec(W)/dt;
      Wprev_vec = B\b;
      Wprev     = reshape(Wprev_vec, I, 2);

      %% B) FIRM
      %----------------------------------------------
      BB = (sig + r - delta + 1/dt) * speye(I) - A(1:I, 1:I);
      bb = profit + J/dt;
      Jprev = BB\bb;

      %% C) NASH BARGAINING
      %----------------------------------------------
      cstar = c(:, 1);
      if gam == 1
         ustar = log(cstar);
      else
         ustar = cstar.^(1-gam)./(1-gam);
      end
      daW   = cstar.^(-gam);
      daJ   = [diff(J)/da; 0];
      dtW   = (W - Wprev)/dt;
      dtJ   = (J - Jprev)/dt;

      if strcmpi(bargain, 'nash')
         zeta1 = ( k^alfa - r*k + daJ.*((r-delta)*a - cstar) + dtJ ) ./ (1-daJ);
         zeta2 = (r-delta)*a - cstar + (ustar - rho*W(:, 2) + dtW(:,1)) ./ daW;
         wprev = bet*zeta1 - (1-bet)*zeta2;
      elseif strcmpi(bargain, 'egalitarian')
         scal = (1-bet)/(rho+sig).*daW + bet/(r-delta+sig).*(1-daJ);
         zeta1 = ( z*k^alfa - r*k + daJ.*((r-delta)*a - cstar) + dtJ ) ./ (r-delta+sig);
         zeta2 = ( ustar + daW.*((r-delta)*a - cstar) - rho*W(:, 2) + dtW(:, 1) ) ./ (rho+sig);
         wprev = (bet*zeta1 - (1-bet)*zeta2) ./ scal;
      end

      %% Save some stuff
      c_t{n}      = c;
      adot_t{n}   = adot;
      A_t{n}      = A;
      w_t{n}      = w;
      W_t{n}      = W;
      J_t{n}      = J;
      profit_t{n} = profit;

      % update values and wage
      J = Jprev;
      W = Wprev;
      w = wprev;
   end

   %                   3.  Iterate Kolmogorov forward
   %---------------------------------------------------------------------
   % Always start from initial distribution revalued by equity price change
   % Trick is: we need integration-preserving interpolation
   anew     = a * (1 + (p0 - steady.p) / (steady.K + steady.p)); % equity gain
   gg_pchip = interp1(anew, steady.gg, a, 'pchip');
   ggreval  = gg_pchip / sum(sum(gg_pchip)) / da;

   % Start Kolmogorov forwars
   gg_t{1}  = vec(ggreval);
   for n = 1:N-1
      AT        = transpose(A_t{n} + Lambda_t{n});
      gg_t{n+1} = (speye(I*2) - dtvec(n)*AT) \ gg_t{n}; % implicit method
   end

   %                      4.  Check market clearing
   %---------------------------------------------------------------------
   p_t(end) = steady.p;
   for n = N:-1:1
      gg_temp = reshape(gg_t{n}, I, 2);

      % unempoyment moves smoothly, vacancies jump
      u_t(n) = sum(gg_temp(:, 2)) * da;
      v_t(n) = theta_t(n) .* u_t(n);

      % Free entry
      FE_t(n) = -xi + q_t(n) * transpose(J_t{n}) * gg_temp(:, 2)/u_t(n) * da;

      % Asset market
      K_t(n)  = (1 - u_t(n)) .* k_t(n);
      d_t(n)  = transpose(profit_t{n}) * gg_temp(:, 1) * da - xi*v_t(n);
      AD_t(n) = a' * sum(gg_temp, 2) * da - K_t(n) - p_t(n); % excess demand for assets

      % previous equity price
      if n > 1
         p_t(n-1) = ( d_t(n-1) + p_t(n)/dtvec(n-1) ) / (1/dtvec(n-1) + r_t(n-1) - delta);
      end
   end

   % Equity price on impact
   EP = p_t(1) - p0;
   toc

   if mod(it, 10) == 0 && it < 100
      dk = 1.5 * dk;
      dtheta = 1.5 * dtheta;
   elseif mod(it, 10) == 0 && it > 100 && it < 150
      dk = 1.2 * dk;
      dtheta = 1.2 * dtheta;
   end

   figure(1)
   subplot(2,2,1); plot(time, FE_t);    hold on; plot(time, zeros(length(time), 1), 'r--');             title('Free entry'); hold off;
   subplot(2,2,2); plot(time, theta_t); hold on; plot(time, steady.theta*ones(length(time), 1), 'r--'); title('tightness');  hold off;
   subplot(2,2,3); plot(time, AD_t);    hold on; plot(time, zeros(length(time), 1), 'r--');             title('Assets');     hold off;
   subplot(2,2,4); plot(time, k_t);     hold on; plot(time, steady.k*ones(length(time), 1), 'r--');     title('capital');    hold off;
   pause(0.001)

   %%                            4. Updating
   %----------------------------------------------------------------------

   % Adjust tightness (increase if firms are too profitable)
   theta_t = theta_t + dtheta .* FE_t;
   disp(['Error in free entry: ', num2str(max(abs(FE_t)))])

   % Adjust interest rate (decrease if HHs demand too much assets)
   k_t = k_t + dk .* AD_t;
   disp(['Error asset market: ', num2str(max(abs(AD_t)))])

   % Adjust equity price on impact (realign with backward iteration)
   p0 = dp0 * p_t(1) + (1-dp0) * p0;
   disp(['Error equity price on impact: ', num2str(EP)])
   disp('----------------------------------------------')

   if max(abs(FE_t)) < tol_ent && max(abs(AD_t(2:end))) < tol_mkt && abs(EP) < tol_eqty
      disp('Found the equilibrium!')
      break
   end

   save('guess_k.mat', 'theta_t', 'k_t', 'p0')
end

%==========================================================================
%                                 Results
%==========================================================================
% Value functions
sol.W_t      = W_t;      % worker value function
sol.J_t      = J_t;      % filled job value function
sol.gg_t     = gg_t;     % worker distribution

% Policy functions
sol.c_t      = c_t;      % consumption
sol.adot_t   = adot_t;   % savings

% Asset market
sol.k_t      = k_t;      % capital per worker
sol.K_t      = K_t;      % aggregate capital
sol.p_t      = p_t;      % equity price
sol.d_t      = d_t;      % dividend
sol.r_t      = r_t;      % dividend
sol.A_t      = A_t;      % transition matrix

% Labor market
sol.u_t      = u_t;      % unemployment rate
sol.v_t      = v_t;      % vacancies
sol.w_t      = w_t;      % wage schedule
sol.f_t      = f_t;      % job-finding rate
sol.q_t      = q_t;      % job-filling rate
sol.theta_t  = theta_t;  % tightness
sol.Lambda_t = Lambda_t; % labor market transition matrix

% TFP
sol.z_t      = z_t;

end
