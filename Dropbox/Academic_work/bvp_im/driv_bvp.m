% this sample driver computes the solution of x' = -x on [0,1] with bc g(x0,xend) =
% x0+xend-1 using the midpoint rule
d=2;
t0=0; tf = pi/2; x0 = [0 0]';; xend = [2 0]';; 
M=100; par=1; prob = 1;
xinit = zeros(M,d);
for n = 1:M
  xinit(n,1) = 0; % 2*sin((n-1)*(tf-t0)/(M-1));
  xinit(n,2) = 0; %2*cos((n-1)*(tf-t0)/(M-1));
end
% use 1 for midpt rule
solver.disc = 1; 
% M is number of solution points
solver.numberofpts = M;
% tolance and iteration max for Newton solver
solver.newtontol = 1E-13;
solver.newtonct = 10;
[X,T,stats] = bvpsol(xinit,t0,tf,x0,xend,d,prob,par,solver);

exsol = zeros(M,2);  exsoldiff = zeros(M,2);
for n = 1:M
  exsoldiff(n,1) = 2*sin((n-1)*(tf-t0)/(M-1))- X(n,1);
  exsoldiff(n,2) = 2*cos((n-1)*(tf-t0)/(M-1))- X(n,2);
end

% norm of exsoldiff is the absolute solution error
