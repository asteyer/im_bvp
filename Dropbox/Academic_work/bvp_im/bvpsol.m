%  function bvpsol implements a 2-point bvp solve using a discretization specified by solver
%  2-point bvp of the form x' = f(x,t), g(x0,xend) = 0
%
%  inputs: xinit, t0, tf, x0, xend, M,d,prob, par, solver, newtontol
%
%  xinit - initial guess for solution ( M x d in size)
%  t0    - left end-point
%  tf    - right end-point
%  x0    - value of left end-point  x(t0) = x0
%  xend  - value of the right end-point x(end) = xend  ; g(x0,xend) = 0
%  d     - dimension of the problem
%  prob  - what problem is being solved
%  par   - vector of problem parameters
%  solverpar - solver parameters; solver.disc is the discretization specified, solver.newtontol = tolerance of newton solver
%
% outputs: X,T,stats
%
% X     - M x d solution vector
% T     - solution time vectors
% stats - various error statistics, stats.newton is the Newton errors 
%
%  discretization of form [J,F] = midptbvpsetup(X,T,H,d,prob,par,M); F = function to be solved, J is its jacobian
function [X,T,stats] = bvpsol(xinit,t0,tf,x0,xend,d,prob,par,solverpar)
  X=xinit;
  M = solverpar.numberofpts;
  T=zeros(M,1); h = (tf-t0)/(M-1);  T(1) = t0; H=zeros(M-1);
  for n=2:M
    T(n) = T(n-1)+h;
    H(n-1) = h;
  end
  if solverpar.disc == 1
    [J,F] = midptbvpsetup(X,T,H,d,prob,par,M);
  end
  Jnorm = norm(J,inf); Fnorm = norm(F,inf);
  if Fnorm > Jnorm
    Jnorm = Fnorm;
  end
  error = norm(F,2)/(1+Jnorm);
  counter = 0; dx = zeros(M,2);
  while error > solverpar.newtontol && counter < solverpar.newtonct 
    dxtemp = -J\F;
    % this loop is not efficient, figure out how to reshape dx correctly;
    for j=1:M
      dx(j,1) = dxtemp(2*(j-1)+1);
      dx(j,2) = dxtemp(2*(j-1)+2);
    end
    X = X+dx;
    [J,F] = midptbvpsetup(X,T,H,d,prob,par,M);
    error = norm(F,2)/(1+Jnorm);
    counter = counter + 1;
  end
  stats.newton = {'error',error,'iterates',counter};
end 
  

