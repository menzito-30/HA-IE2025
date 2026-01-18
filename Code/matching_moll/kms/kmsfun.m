function st = kmsfun(param, gridd, num, guess, bargain)
% Function to calculate the stationary equilibrium of a Krusell-Mukoyama-Sahin model.
%
% Input:
%     param : struct, contains structural parameters of the model
%     gridd : struct, specifies grid for assets
%     num   : struct, contains numerical parameters
%     guess : struct, initial guess for value functions, capital, interest rate. wage schedule
%
% Output:
%     st   : struct, model objects in stationary equilibrium

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

I         = gridd.I;
amin      = gridd.amin;
amax      = gridd.amax;

tol_hjb   = num.tol_hjb;
tol_ent   = num.tol_ent;
tol_mkt   = num.tol_mkt;
maxit_hjb = num.maxit_hjb;
maxit_out = num.maxit_out;
Delta     = num.Delta;
dtheta    = num.dtheta;
dk        = num.dk;

theta     = guess.theta;
k         = guess.k;

% Finish building the grid
a    = linspace(amin, amax, I)';
aa   = [a, a]; % to work in old MATLAB
da   = (amax-amin) / (I-1);

% Finish initial guess
r       = alfa * k ^ (alfa-1);
w       = ones(I, 1) * bet * (k^alfa - r*k);
W       = zeros(I, 2);
if gam == 1 % log utility
   V(:, 1) = log(w + (r-delta)*a) / rho; % keep assets constant
   V(:, 2) = log(h + (r-delta)*a) / rho;
else % CRRA
   V(:, 1) = (w + (r-delta)*a) .^ (1-gam) / (1-gam) / rho;
   V(:, 2) = (h + (r-delta)*a) .^ (1-gam) / (1-gam) / rho;
end
J       = (k^alfa - r*k - w) / (sig + r - delta); % disregard change in assets

% Preallocation
WF = zeros(I, 2);
WB = zeros(I, 2);
Ak = cell(2, 1);
FE = zeros(maxit_out, 1);
AD = zeros(maxit_out, 1);

%==========================================================================
%                                Iteration
%==========================================================================
for it = 1:maxit_out
   disp(['Outer loop, Iteration ' int2str(it)]);

   % update functions of tightness
   f = chi * theta^(1-eta);
   q = chi * theta^(-eta);
   u = sig / (sig + f);
   v = theta * u;

   % labor market transition matrix
   lam    = [-sig, sig; f, -f];
   Lambda = [speye(I)*lam(1,1), speye(I)*lam(1,2);
             speye(I)*lam(2,1), speye(I)*lam(2,2)];

   % update functions of r and w
   r = alfa * k ^ (alfa-1);     % interest rate
   K = (1-u) * k;               % aggregate capital
   y = [w, ones(I, 1)*h];       % income
   profit = k^alfa - r*k - w;

   % Not updating value function on purpose. Better then initial guess.

   %                    1. Solve worker HJB equation
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
      c    = If.*cF + Ib.*cB + I0.*c0;
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
      B1 = (rho + 1/Delta)*speye(I*2) - A - Lambda;
      b1 = vec(util) + vec(W)/Delta;
      Wnew_vec = B1\b1;
      Wnew     = reshape(Wnew_vec, I, 2);

      % check convergence
      con = max(max(abs(Wnew - W)));
      if con < tol_hjb
         disp(['Value Function W Converged, Iteration = ', num2str(i)]);
         break
      end

      % update
      W = Wnew;
   end

   %                 2. Solve Kolmogorov Forward Equation
   %----------------------------------------------------------------------
   AT  = (A + Lambda)';
   nul = zeros(I*2, 1);

   % Normalize g(i,s) in one point (ensure later that it integrates to 1)
   i_fix        = 1;
   nul(i_fix)   = 1;
   row          = [zeros(1, i_fix-1), 1, zeros(1, 2*I-i_fix)];
   AT(i_fix, :) = row;

   % Solve linear system
   gvec = AT\nul;
   gsum = sum(gvec) * da;
   g    = gvec./gsum;
   gg   = reshape(g, I, 2);


   %                     3. Solve firm HJB equation
   %----------------------------------------------------------------------
   % LoM for savings, exogenous for firm, can transition matrix M from above
   B2 = (1/Delta + sig + r - delta)*speye(I) - A(1:I, 1:I);
   for i = 1:maxit_hjb
      % set up and solve linear system
      b2 = profit + J/Delta;
      Jnew = B2\b2;

      % check convergence
      con = max(max(abs(Jnew - J)));
      if con < tol_hjb
         disp(['Value Function J Converged, Iteration = ', num2str(i)]);
         break
      end

      % update
      J = Jnew;
   end

   %                 4. Market clearing and wage setting
   %----------------------------------------------------------------------
   % Free entry
   FE(it) = -xi + q * J' * gg(:, 2)/u * da;

   % Asset mkt
   d      = profit' * gg(:, 1) * da - xi*v;
   p      = d / (r-delta);
   AD(it) = a' * sum(gg, 2) * da - K - p; % excess demand of assets

  % Bargaining
   cstar = c(:, 1);
   if gam == 1
      ustar = log(cstar);
   else
      ustar = cstar.^(1-gam)./(1-gam);
   end
   dW = cstar.^(-gam);
   dJ = [diff(J)/da; 0];
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

   %                             4. Updating
   %----------------------------------------------------------------------
   % make updating more aggresive over time
   if mod(it, 50) == 49
      dtheta = dtheta * 1.3;
   end

   % Adjust tightness (increase if firms are too profitable)
   theta = theta + dtheta * FE(it);
   disp(['Error in free entry: ', num2str(FE(it))])
   disp(['New tightness: ', num2str(theta)])

   % Adjust interest rate (decrease if HHs demand too much assets)
   k = k + dk * AD(it);
   disp(['Error asset market: ', num2str(AD(it))])
   disp(['New capital: ', num2str(k)])
   disp('----------------------------------------------')

   % Check if we're in equilibrium
   if abs(FE(it)) < tol_ent && abs(AD(it)) < tol_mkt
      disp(['Found the equilibrium! theta = ', num2str(theta), 'r = ', num2str(r)])
      break
   end
end

%==========================================================================
%                                 Results
%==========================================================================
% Value functions
st.W      = W;      % worker value function
st.J      = J;      % filled job value function
st.c      = c;      % consumption policy
st.adot   = adot;   % saving policy
st.gg     = gg;     % worker distribution

% Asset market
st.k      = k;      % capital per worker
st.K      = K;      % aggregate capital
st.p      = p;      % equity price
st.d      = d;      % dividend
st.r      = r;      % interest rate

% Labor market
st.u      = u;      % unemployment rate
st.v      = v;      % vacancies
st.w      = w;      % wage schedule
st.f      = f;      % job-finding rate
st.q      = q;      % job-filling rate
st.theta  = theta;  % tightness
st.Lambda = Lambda; % labor market transition matrix

% Statistics
st.k_return  = 1 - exp(-4*(r-delta));     % annualized return on capital
st.i_share   = delta * k^(1-alfa);        % investment share
st.meanasset = a' * sum(gg, 2) * da;      % mean wealth
st.meanwage  = w' * gg(:,1) * da / (1-u); % mean wage
st.minwage   = min(w);                    % minimum wage
st.meanmin   = st.meanwage / st.minwage;  % mean-min ratio
st.repratio  = h / st.meanwage;           % mean replacement ratio
st.output    = (1-u) * k^alfa;            % output

% mean surplus from a match
st.surplus = (W(:, 1) - W(:, 2) + J)' * gg(:, 1) * da;

end