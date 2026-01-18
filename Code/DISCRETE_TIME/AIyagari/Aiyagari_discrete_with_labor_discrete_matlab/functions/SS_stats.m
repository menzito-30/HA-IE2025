%% Calculating SS statistics
% --------------------------

SS.K = grid.a * sum(JD,2);
SS.C = sum(sum(pol.c .* JD));
SS.N = sum(sum(mesh.z .* pol.l .* JD));
SS.Y = SS.K^par.alpha * SS.N^(1-par.alpha);
SS.ProbBoC = sum(JD(1,:));    % Mass at Borrowing constraint grid.a(1)
SS.Borrow = sum(JD(mesh.a<0));  % mass with negative savings
SS.Save = sum(JD(mesh.a>0));    % mass with positive savings
SS.KN = SS.K/SS.N;