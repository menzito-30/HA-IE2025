%==========================================================================
%         Calibratring Vacancy Cost Krusell-Mukoyama-Sahin model
%
%
% Author: Bence Bardoczy
% Date  : 5/31/2017
%==========================================================================

%=========================================================================
%%                           Housekeeping
%=========================================================================
close all;
clear all;
clc;

%=========================================================================
%%                         Initialization
%=========================================================================
% Structural parameters, unit of time is quarter
gam   = 1;      % relative risk aversion
rho   = 0.01;   % discount rate
chi   = 1.7935; % matching efficiency
sig   = 0.1038; % separation rate
eta   = 0.72;   % matching elasticity to unemp.
bet   = 0.72;   % worker bargaining power
h     = 0.75;   % home production
alfa  = 0.3;    % capital share
delta = 0.021;  % depreciation rate
theta = 1;      % tightness

% Numerical parameters
tol_hjb   = 1e-8; % tolerance level for HJB equations
tol_ent   = 1e-6; % tolerance level for free entry condition
tol_mkt   = 1e-6; % tolerance level for asset market clearing
maxit_hjb = 50;   % maximum number of iterations when solving HJB equation
maxit_out = 2000; % maximum number of iterations in outer loop
Delta     = 2000; % step size in HJB equation
dxi       = 1e-1; % step size in updating vacancy cost
dk        = 2e-3; % step size in updating interest rate

% Bargaining solution
% bargain = 'nash';
bargain = 'egalitarian'

% Asset grid
I    = 1001;
amin = 0;
amax = 150;
a    = linspace(amin, amax, I)';
aa   = [a, a]; % to work in old MATLAB
da   = (amax-amin) / (I-1);

% Initial guesses
xi = 0.132;
k  = 1.2 * (alfa / (rho+delta)) ^ (1/(1-alfa));
r  = alfa * k^(alfa-1);
w  = ones(I, 1) * bet * (k^alfa - r*k);

% Implied variables
f = chi * theta^(1-eta);
q = chi * theta^(-eta);
u = sig / (sig + f);
v = theta * u;

% labor market transition matrix
lam    = [-sig, sig; f, -f];
Lambda = [speye(I)*lam(1,1), speye(I)*lam(1,2);
          speye(I)*lam(2,1), speye(I)*lam(2,2)];

% Initial guesses for HJB equations
W = zeros(I, 2);
if gam == 1 % log utility
   W(:, 1) = log(w + (r-delta)*a) / rho; % keep assets constant
   W(:, 2) = log(h + (r-delta)*a) / rho;
else % CRRA
   W(:, 1) = (w + (r-delta)*a) .^ (1-gam) / (1-gam) / rho;
   W(:, 2) = (h + (r-delta)*a) .^ (1-gam) / (1-gam) / rho;
end
W0 = W;
J = (k^alfa - r*k - w) / (sig + r - delta); % disregard change in assets

% Preallocation
VF = zeros(I, 2);
VB = zeros(I, 2);
Ak = cell(2, 1);
FE = zeros(maxit_out, 1);
AD = zeros(maxit_out, 1);

%=========================================================================
%%                            Iteration
%=========================================================================
tic
for it = 1:maxit_out
   disp(['Outer loop, Iteration ' int2str(it)]);

   % update functions of k and w
   r = alfa * k ^ (alfa-1);     % interest rate
   K = (1-u) * k;               % aggregate capital
   y = [w, ones(I, 1)*h];       % income
   profit = k^alfa - r*k - w;

   % Not updating value function on purpose. Better then initial guess.
   % V = V0;

   %%                  1. Solve worker HJB equation
   %---------------------------------------------------------------------
   for i = 1:maxit_hjb
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

      % diagonals
      X = -min(sB, 0) / da;
      Z =  max(sF, 0) / da;
      Y = -X - Z;

      % prepare for spdiags
      X = X(2:I,:);
      Z = [zeros(1,2); Z(1:I-1,:)];

      % transition matrix
      for j = 1:2
         Ak{j} = spdiags(Y(:,j), 0, I, I) + spdiags(X(:,j), -1, I, I) + spdiags(Z(:,j), 1, I, I);
      end
      M = [Ak{1}, sparse(I,I); sparse(I,I), Ak{2}];
      A = M + Lambda;

      % set up and solve linear system
      B = (rho + 1/Delta)*speye(I*2) - A;
      b = vec(util) + vec(W)/Delta;
      Wnew_vec = B\b;
      Wnew = reshape(Wnew_vec, I, 2);

      % check convergence
      Wchange = Wnew - W;
      con = max(max(abs(Wchange)));
      if con < tol_hjb
         disp(['Value Function V Converged, Iteration = ', num2str(i)]);
         break
      end

      % update
      W = Wnew;
   end

   %%                 2. Solve Kolmogorov Forward Equation
   %----------------------------------------------------------------------
   AT  = A';
   nul = zeros(I*2, 1);

   % Normalize g(i,j) in one point (ensure later that it integrates to 1)
   i_fix        = 1;
   nul(i_fix)   = 1;
   row          = [zeros(1, i_fix-1), 1, zeros(1, 2*I-i_fix)];
   AT(i_fix, :) = row;

   % Solve linear system
   gvec = AT\nul;
   gsum = sum(gvec) * da;
   g    = gvec./gsum;
   gg   = reshape(g, I, 2);


   %%                    3. Solve firm HJB equation
   %----------------------------------------------------------------------
   % LoM for savings, exogenous for firm, can transition matrix M from above
   BB = (1/Delta + sig + r - delta)*speye(I) - M(1:I, 1:I);
   for i = 1:maxit_hjb
      % set up and solve linear system
      bb = vec(profit) + vec(J)/Delta;
      Jnew = BB\bb;

      % check convergence
      Jchange = Jnew - J;
      con = max(max(abs(Jchange)));
      if con < tol_hjb
         disp(['Value Function J Converged, Iteration = ', num2str(i)]);
         break
      end

      % update
      J = Jnew;
   end

   %%                 4. Market clearing and wage setting
   %----------------------------------------------------------------------
   % Free entry
   FE(it) = -xi + q * J' * gg(:, 2)/u * da;

   % Asset mkt
   div    = profit' * gg(:, 1) * da - xi*v;
   p      = div / (r-delta);
   AD(it) = a' * sum(gg, 2) * da - K - p; % excess demand of assets

  % Nash bargaining
   cstar = c(:, 1);
   if gam == 1
      ustar = log(cstar);
   else
      ustar = cstar.^(1-gam)./(1-gam);
   end
   dW    = cstar.^(-gam);
   dJ    = [diff(J)/da; 0];
   if strcmpi(bargain, 'nash')
      zeta1 = ( k^alfa - r*k + dJ.*((r-delta)*a - cstar) ) ./ (1-dJ);
      zeta2 = (r-delta)*a - cstar + (ustar - rho*W(:, 2)) ./ dW;
      w     = bet*zeta1 - (1-bet)*zeta2;
   elseif strcmpi(bargain, 'egalitarian')
      scal = (1-bet)/(rho+sig).*dW + bet/(r-delta+sig).*(1-dJ);
      zeta1 = ( k^alfa - r*k + dJ.*((r-delta)*a - cstar) ) ./ (r-delta+sig);
      zeta2 = ( ustar + dW.*((r-delta)*a - cstar) - rho*W(:, 2) ) ./ (rho+sig);
      w     = (bet*zeta1 - (1-bet)*zeta2) ./ scal;
   end

   %%                            4. Updating
   %----------------------------------------------------------------------
   % Adjust vacancy cost (increase if firms are too profitable)
   if abs(FE(it)) > tol_ent
      xi = xi + dxi * FE(it);
      disp(['Error in free entry: ', num2str(FE(it))])
      disp(['New vacancy cost: ', num2str(xi)])
   end
   % Adjust interest rate (decrease if HHs demand too much assets)
   if abs(AD(it)) > tol_mkt
      k = k + dk * AD(it);
      disp(['Error asset market: ', num2str(AD(it))])
      disp(['New interest rate: ', num2str(r)])
   end
   disp('----------------------------------------------')

   if abs(FE(it)) < tol_ent && abs(AD(it)) < tol_mkt
      disp(['Found the equilibrium! xi = ', num2str(xi), 'r = ', num2str(r)])
      break
   end
end
toc
