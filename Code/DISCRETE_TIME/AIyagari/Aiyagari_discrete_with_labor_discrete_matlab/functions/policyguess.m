function c_guess = policyguess(mpar, mesh, R, W)

l0 = 0.6 * ones(mpar.na, mpar.nz);  % initial guess: l=0.6 in all states
c_guess = W*mesh.z.*l0 + (R-1)*mesh.a;  % initial guess: a'=a in all states
assert(all(c_guess(:)>0),'Error policy guess: not all c_guess>0');
