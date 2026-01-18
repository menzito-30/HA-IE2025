
% grid for idiosyncratic income process
% ---------------------------------------
[grid.z, Pz] = rouwen(par.rhoz, par.uncond_muz, par.uncond_sigmaz, mpar.nz);
grid.z = exp(grid.z);    % take levels

% ergodic distribution
% ---------------------------------------
grid.invdistr = Pz^1000;
grid.invdistr = grid.invdistr(1,:);

% guess aggregate capital-to-labor ratio
% ---------------------------------------
grid.N = grid.invdistr(1,:) * grid.z;    % aggregate Employment if everybody works fully
grid.K = ( (par.alpha * grid.N^(1-par.alpha)) / (1/par.beta - 1 + par.delta) )^(1/(1-par.alpha)); % RepAgent SS
grid.KN = grid.K/grid.N;

% grid for individual capital
% ---------------------------------------
grid.a_min = -1.5;
grid.a_max = 10*grid.K;
grid.a = (exp(exp(linspace(0,log(log(grid.a_max - grid.a_min+1)+1),mpar.na))-1)-1+grid.a_min);

% state space (a,z)
% ---------------------------------------
[mesh.a, mesh.z] = ndgrid(grid.a,grid.z);