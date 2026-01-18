function pol = policies_EGM(par, mpar, grid, mesh, c_guess, Pz, R, W)

distEGM = 9999;
countEGM = 0;

while distEGM>mpar.crit && countEGM<10000
    countEGM=countEGM+1;

    % compute marginal utility of consumption (t+1)
	% ---------------------------------------------
    mu_c = c_guess.^(-par.sigma);

    % compute right hand side of Euler equation
	% ---------------------------------------------
    rhs_aux = par.beta * R * mu_c * Pz';

    % compute candidate consumption c(t)
	% ---------------------------------------------
    c_cand = rhs_aux.^(-1/par.sigma);   % c_cand = c(a',z)
    assert(all(c_cand(:)>0),'Error EGM: not all c_cand>0');
    
    % compute candidate labor l(t)
	% ---------------------------------------------
    l_cand = (1/par.psi1 * W * mesh.z .* c_cand.^(-par.sigma)).^(1/par.psi2);
    l_cand(l_cand<0)=0; l_cand(l_cand>1)=1;
    
    % compute candidate savings policy a'(t)
	% ---------------------------------------------
    a_cand = (c_cand + mesh.a - W.*mesh.z.*l_cand) ./ R;   % a_cand=a(t)
    
    % compute threshold a, that leads to a'(t)=a_min (given z)
	% ---------------------------------------------
    a_thresh = a_cand(1,:); 
    a_thresh = repmat(a_thresh, [mpar.na 1]);
    % check in which states a (for given z) constraint binds
    Constrained = ( mesh.a < a_thresh );

    % interpolate candidate policy functions on grid
	% ---------------------------------------------
	c_pol=NaN(mpar.na, mpar.nz);
	l_pol=NaN(mpar.na, mpar.nz);
	a_pol=NaN(mpar.na, mpar.nz);
    for iz=1:mpar.nz
		c_pol(:,iz) = interp1(a_cand(:,iz),c_cand(:,iz),grid.a,'linear','extrap');
        a_pol(:,iz) = interp1(a_cand(:,iz),grid.a      ,grid.a,'linear','extrap');
        l_pol(:,iz) = interp1(a_cand(:,iz),l_cand(:,iz),grid.a,'linear','extrap');
		
		    
        % correct for constrained states
		% ---------------------------------------------
		izConstrained = Constrained(:,iz);
        if any(izConstrained)
            [c_pol(izConstrained,iz), l_pol(izConstrained,iz)] = egm_solve_constrained(par, grid, grid.a(izConstrained), iz, R, W);
            a_pol(izConstrained,iz) = grid.a(1);
        end
    end

    % check convergence
	% ---------------------------------------------
    distEGM = max(abs( c_pol(:)-c_guess(:) ));
    % update guess
	% ---------------------------------------------
    c_guess = c_pol;
    
end

% assign output
% ---------------------------------------
pol.c = c_pol;
pol.a = a_pol;
pol.l = l_pol;
pol.constrained=Constrained;

end



% _________________________________________________________________________
%% EGM solver when asset choice is constrained
% _________________________________________________________________________
function [c,l] = egm_solve_constrained(par,grid,a,iz,R,W)

% check whether a is a column vector, if not make it one
% ---------------------------------------
if ~iscolumn(a); a=a'; end

% initial guess for l and c
% ---------------------------------------
l = 0.6 * ones(size(a));
c = R*a + W*grid.z(iz)*l - grid.a(1);

% Rootfinding for l via Newton
% ---------------------------------------
critEGM = all(abs(W*grid.z(iz)*c.^(-par.sigma) - par.psi1 * l.^par.psi2));
countEGM=0; tolEGM=1e-7;
while critEGM>tolEGM && countEGM<1000
    countEGM = countEGM + 1;
    
	% intra-period FOC for labor supply (function value)
	% ---------------------------------------
    f = W*grid.z(iz)*c.^(-par.sigma) - par.psi1 * l.^par.psi2;

    % derivative of intra-period FOC for labor supply
	% ---------------------------------------
    J = -par.sigma*(W*grid.z(iz))^2 *c.^(-par.sigma-1) - par.psi1*par.psi2 * l.^(par.psi2-1);
	
	% update labor supply choice
	% ---------------------------------------
    lnew = l - f./J;
    % update consumption
	% ---------------------------------------
    cnew = R*a + W*grid.z(iz)*l - grid.a(1);
    
    % compute critical value
	% ---------------------------------------
    critEGM = max(abs( lnew(:)-l(:) ));
    l=lnew; c=cnew;
    
end

% check validity of solution
% ---------------------------------------
assert(countEGM<1000, 'Error EGM: egm_solve_constrained did not converge')

% impose l in (0,1)
% ---------------------------------------
l(l<0)=0; l(l>1)=1;
% compute consumption policy
% ---------------------------------------
c = R*a + W*grid.z(iz)*l - grid.a(1);

end
