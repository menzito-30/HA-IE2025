function [JD, countJD, distJD] = JDiteration(mpar, grid, Pz, pol)
% _________________________________________________________________________
%% Construct Markov transition matrix
% _________________________________________________________________________
% Find next-smallest on-grid value (and associated weight) for pol.a(a,z)
% -------------------------------------------------
[dista, ida] = genweight(pol.a, grid.a);        % get w_i(a,z) and i(a,z)

% state in (t) ("where we come from")
% -------------------------------------------------
rowindex = repmat(1:mpar.na*mpar.nz, [1 2*mpar.nz]);

% state in (t+1) ("where we go to")
% -------------------------------------------------
ida = repmat(ida(:),[1 mpar.nz]);   % get index-number for state (a,z)
idz=kron(1:mpar.nz, ones(1,mpar.na*mpar.nz));
colindex1 = sub2ind([mpar.na mpar.nz], ida(:), idz(:));    % get column index for a'(a,z)=a_i
colindex2 = sub2ind([mpar.na mpar.nz], ida(:)+1, idz(:));  % get column index for a'(a,z)=a_i+1

% probability of transition ("prob of going there")
% -------------------------------------------------
Hz = kron(Pz, ones(mpar.na,1));
weight1 = (1-dista(:)) .* Hz;
weight2 = dista(:) .* Hz;

% sparse transition matrix
% -------------------------------------------------
H = sparse(rowindex, [colindex1 colindex2], [weight1 weight2], ...
        mpar.na*mpar.nz, mpar.na*mpar.nz);
    
% _________________________________________________________________________
%% Compute stationary distribution via iteration
% _________________________________________________________________________
% get initial guess
% -------------------------------------------------
[JD,~] = eigs(H',1);  % JD = eigenvector to H associated w/ largest eigenvalue
JD=JD'./sum(JD);

% Compute Stationary Distribution
% -------------------------------------------------
distJD=9999;
countJD=1;
while (distJD>1e-14 || countJD<50) && countJD<10000 %at least 50 iterations
    countJD=countJD+1;

    % iterate one period/step
	% -----------------------
    JD_next = JD * H;
    JD_next = full(JD_next);
    JD_next = JD_next./sum(JD_next);
    
	% check convergence
	% -----------------------
    distJD=max((abs(JD_next(:)-JD(:))));
    % update distribution
	% -----------------------
    JD=JD_next;
end

% reshape to (a,z)
% -------------------------------------------------
JD = reshape(JD,[mpar.na mpar.nz]);

% check validity of solution
% -------------------------------------------------
assert(countJD<10000, 'Error JD: JDiteration did not converge'); 

end



% _________________________________________________________________________
%% Finding the nearest neighbor for a given policy rule
% _________________________________________________________________________
function [ weight,index ] = genweight( x,xgrid )
% function: GENWEIGHT generates weights and indexes used for linear interpolation
%
% X: Points at which function is to be interpolated.
% xgrid: grid points at which function is measured
% no extrapolation allowed
[~,index] = histc(x,xgrid);
index(x<=xgrid(1))=1;
index(x>=xgrid(end))=length(xgrid)-1;

weight = (x-xgrid(index))./(xgrid(index+1)-xgrid(index)); % weight xm of higher gridpoint
weight(weight<=0) = 1.e-16; % no extrapolation
weight(weight>=1) = 1-1.e-16; % no extrapolation

end