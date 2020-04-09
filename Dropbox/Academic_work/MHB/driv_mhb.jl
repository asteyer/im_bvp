using Pkg
using PyPlot
#############################################
# Basic driver script to run the MHB code 
#############################################
#########################################################
# Example driver file to approximate the periodic orbit 
# from page 162 of Hirsch, Smale, Devaney
#  unique hyperbolic periodic orbit:
#   (x,y) = (cos(t),sin(t))
#########################################################

using LinearAlgebra

### Must include the following files
include("mhb.jl");
include("ode.jl");
include("quadrature.jl");

### Set number of harmonics ###
M    = 1;
## Set the exact period ################ 
T    = 2 * pi;   
#### set dimension ###
dim  = 2;

####### set the odepars
problem = 2;  ### Problem = 2 is the nonlinear problem, you can add your own problem as desired.
modelpars = [];
odepar = ode.odepar(dim,problem,modelpars);

#### set the solver pars
tol     = 1e-15;
maxiter = 10;
### set number of points in Fourier coefficient integral
numpts  = 200;
### Solve for the period
varyT   = true;

### First test, Test 3.1, initialize with true solution yields zero iterations 
### and machine precision err ###
### Initialize with exact solution (x,y) = (cos(t),sin(t))
A = zeros(dim,2*M+1)
A[:,3] = [.32; .73];
A[:,2] = [.88 ; .11];
A[:,1] = [.1 ; -.1];
T      = 6;
println("Starting Newton iteration")
(A,T,err,iter) = mhb.solveResidualNewton(A,M,T,mhbpar);
println("Newton iteration completed")

#### The following code sets up a vector of solutions to plot
tf = 10; ## final time
N = 100; ## number of time points
dt = tf/(N-1); ## 

X    = zeros(dim,N); ## state space solution 
Tvec = zeros(N);     ## time vector corresponding to tf,N,dt
for m = 1:N
  Tvec[m]   = (m-1)*dt;
  ccs    = mhb.getCosSinVec(Tvec[m],T,M);
  p      = A*ccs;
  X[:,m] = p[:]
end

println("fin")
