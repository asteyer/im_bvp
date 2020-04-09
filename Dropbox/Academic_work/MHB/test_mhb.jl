using Pkg
using PyPlot
using LinearAlgebra
include("mhb.jl");
include("ode.jl");

### Choose which tests to run
test1 = true;   ## Test if linear residual is correct
test2 = true;   ## Test MHB solver for 2D linear ODE
test3 = true;   ## Test MHB solver for 2D nonlinear ODE

if test1 ## then test residual for linear ODE u' = [0 w ; -w 0]*u, w > 0
         ## Soutions are of the form: 
         ## u =  a1 * [1 ; 1im]*exp(w*1im*t) + a2 * [1 ; -1im]*exp(-w*1im*t);

  #### M = 1 since only need first harmonic ###
  M    = 1;
  #### Time for the test ####
  time = 3.274;
  w    = 2.3;  # w = 2*pi /T 
  # T  is the period
  T    = 2 * pi /w; 

  ### dim is the dimension of the problem ###
  dim  = 2;

  #### Initialize the Fourier coefficient array
  A = zeros(dim,2*M+1)
  ### solutions are of the form 
  ##  a1 * [-cos(w*t) ; sin(w*t)] + a2 * [sin(w*t) ; cos(w*t)]
  a1 = 1.2; a2 = -3;
  A[:,3] = [a2; a1];
  A[:,2] = [-a1 ; a2];

  #### Set the problem ####
  problem   = 1;
  ### set the model pars ###
  modelpars = [w];
  ### set the vector field pars ###
  odepar = ode.odepar(dim,problem,modelpars);
  #### compute the residual  ####
  R = mhb.getMHBResidual(A,time,M,T,odepar);

  err = norm(R,Inf);
  if err > 1e-12
    print("Test 1 fails, err = "); println(err);
  else
    print("Test 1 passes, err = "); println(err);
  end

end

if test2 ## then test Newton iteration for 2D linear ode

  ##########################################################################
  #### Test convergence of Newton iteration for the linear ode:
  ####               u' = [0 w ; -w 0]*u, w > 0
  #### Solutions are of the form 
  ####      u = a1*[1 ; 1im]*exp(w*1im*t) + a2*[1  ; -1im]*exp(-w*1im*t);
  ##########################################################################

  #### M = 1 since only need first harmonic ###
  M    = 1;

  #### M = 1 since only need first harmonic ###
  M    = 1;
  #### Time for the test ####
  time = 3.274;
  ### Specify the value of w ###
  w    = 2.3;  # w = 2*pi /T 
  # T is the true period
  T    = 2 * pi /w; 

  ### dim is the dimension of the problem ###
  dim  = 2;

  #### Initialize the Fourier coefficient array
  A = zeros(dim,2*M+1)
  ### solutions are of the form 
  ##  a1 * [-cos(w*t) ; sin(w*t)] + a2 * [sin(w*t) ; cos(w*t)]
  a1 = 1.2; a2 = -3;
  A[:,3] = [a2; a1];
  A[:,2] = [-a1 ; a2];

  #### Set the problem ####
  problem   = 1;
  ### set the model pars ###
  modelpars = [w];
  ### set the vector field pars ###
  odepar = ode.odepar(dim,problem,modelpars);
  #### compute the residual  ####
  R = mhb.getMHBResidual(A,time,M,T,odepar);
  ### Set solver tolerance and max no. of iterations
  tol     = 1e-15;
  maxiter = 10;
  ### Set no. of integration points in the inverse Fourier transform
  numpts  = 200;
  ### Solve for the period
  varyT   = true;

  ## Initialize ODE paramters
  problem   = 1;
  modelpars = [w];
  odepar    = ode.odepar(dim,problem,modelpars);
  ### Initialize MHB paramters
  mhbpar         = mhb.mhbpar(odepar,tol,maxiter,numpts,varyT);

  ### First test, Test 2.1, initialize with true solution yields zero iterations 
  ### and machine precision err ###
 (A,T,err,iter) = mhb.solveResidualNewton(A,M,T,mhbpar);
  if (err > 1e-13) || (abs(T-2*pi/w) > 1E-13) || iter > 0
    print("Test 2.1 fails, err, iterations, period err ="); print(err); print(", "); print(iter); print(", "); println(abs(T-2*pi/w));
  else
    print("Test 2.1 passes, err, iterations, period err ="); print(err); print(", "); print(iter); print(", "); println(abs(T-2*pi/w));
  end

  ### Second test, Test 2.2, Newton iteration converges 
  ### when initial guess is near a true solution
  #### 
  #### Set a perturbed initial conditoin
  a1 = 1.2; a2 = -3;
  A[:,3] = [(1-.1)*a2; (1+.2)*a1];
  A[:,2] = [-(1-.3)*a1 ; (1+.09)*a2];
  T      = (1+.3)*T;
 (A,T,err,iter) = mhb.solveResidualNewton(A,M,T,mhbpar);
  if (err > 1e-13) || (abs(T-2*pi/w) > 1E-13) 
    print("Test 2.2 fails, err, iterations, period err ="); print(err); print(", "); print(iter); print(", "); println(abs(T-2*pi/w));
  else
    print("Test 2.2 passes, err, iterations, period err ="); print(err); print(", "); print(iter); print(", "); println(abs(T-2*pi/w));
  end

end

if test3 ## then test Newton iteration for 2D nonlinear ode
         ##  from page 162 of Hirsch, Smale, Devaney
         ##  unique hyperbolic periodic orbit:
         ##   (x,y) = (cos(t),sin(t))
         ###############################################

  ### Set number of harmonics ###
  M    = 1;
  ## Set the exact period ################ 
  T    = 2 * pi;   
  #### set dimension ###
  dim  = 2;


  ####### set the odepars
  problem = 2;
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
  A[:,3] = [0; 1];
  A[:,2] = [1 ; 0];
  A[:,1] = [0 ; 0];
  (A,T,err,iter) = mhb.solveResidualNewton(A,M,T,mhbpar);

  if (err > 1e-13) || (abs(T-2*pi) > 1E-13) || iter > 0
    print("Test 3.1 fails, err, iterations, period err ="); print(err); print(", "); print(iter); print(", "); println(abs(T-2*pi));
  else
    print("Test 3.1 passes, err, iterations, period err ="); print(err); print(", "); print(iter); print(", "); println(abs(T-2*pi));
  end

  ### Second test, Test 3.2, Newton iteration converges 
  ### when initial guess is near a true solution
  #### 
  #### Set a perturbed initial conditoin
  A[:,3] = [-.1; 1+.3];
  A[:,2] = [.88 ; .95];
  A[:,1] = [0.01 ; -.03];
  T      = (1-.25)*2*pi;
 (A,T,err,iter) = mhb.solveResidualNewton(A,M,T,mhbpar);
  if (err > 1e-13) || (abs(T-2*pi) > 1E-13) 
    print("Test 3.2 fails, err, iterations, period err ="); print(err); print(", "); print(iter); print(", "); println(abs(T-2*pi));
  else
    print("Test 3.2 passes, err, iterations, period err ="); print(err); print(", "); print(iter); print(", "); println(abs(T-2*pi));
  end

end


println("fin")
