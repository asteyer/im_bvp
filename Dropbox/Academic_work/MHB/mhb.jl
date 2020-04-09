module mhb

  include("ode.jl")
  using LinearAlgebra
  using QuadGK
  include("quadrature.jl")

  ##############################################################
  # returns a vector of the following form 
  # (a0,a1,...,a_M,b_1,...,b_M)^T 
  # where ansatz is a0 + sum_{j=0}^{M}( a_j * cos(2*pi*j*t/T) + 
  # b_j * sin(2*pi*j*t/T) ) 
  ##############################################################
  function getCosSinVec(time,T,M)
    cosvec = zeros(M,1);
    sinvec = zeros(M,1);
    for j = 1:M
      cosvec[j,1] = cos(2*pi*time*j/T);
      sinvec[j,1] = sin(2*pi*time*j/T);
    end
    return [1; cosvec; sinvec];
  end

  #############################################################
  # returns vector of coefficients from time differentiation
  # of as computed via getCosSinVec:
  #############################################################
  function getCosSinVecDot(time,T,M)
    cosvecdot = zeros(M,1);
    sinvecdot = zeros(M,1);
    for j =1:M
      cosvecdot[j] = -2*pi*j/T;
      sinvecdot[j] =  2*pi*j/T;
    end
    return [0 ; cosvecdot ; sinvecdot];
  end

  ##############################################################
  # This function returns the residual p' - f(p,time,par)
  # where p = a0 + sum_{j=0}^{M}( a_j * cos(2*pi*j*t/T) + 
  # b_j * sin(2*pi*j*t/T) ) 
  ##############################################################
  function getMHBResidual(A,time,M,T,odepar)
    ccs             = getCosSinVec(time,T,M);
    ccsdot          = zeros(2*M+1,1);
    ccsdot[2:M+1]   = ccs[M+2:end];
    ccsdot[M+2:end] = ccs[2:M+1];
    ccsdotcoeff     = getCosSinVecDot(time,T,M);
    ccsdot          = ccsdotcoeff .* ccsdot; 
    p           = A*ccs;
    pdot        = A*ccsdot; 
    return (pdot - ode.frhs(p,time,odepar));
  end 

  ##############################################################
  # This function returns the Fourier coefficients of
  # coefficients of p' - f(p,time,par) 
  # where   # where p = a0 + sum_{j=0}^{M}( a_j * cos(2*pi*j*t/T) + 
  # b_j * sin(2*pi*j*t/T) ) 
  ##############################################################
  function getMHBResidualFC(A,M,T,mhbpar)
    odepar    = mhbpar.odepar;
    dim       = mhbpar.odepar.dim;
    numpts    = mhbpar.numpts;
    Rint      = zeros(dim,2*M+1);
    f0        = y -> mhb.getMHBResidual(A,y,M,T,odepar);
    Rint[:,1] = quadrature.integrate(0,T,numpts,f0,dim);
    for m = 1:M
      fcos          = y -> mhb.getMHBResidual(A,y,M,T,odepar)*cos(2*pi*y*m/T);
      fsin          = y -> mhb.getMHBResidual(A,y,M,T,odepar)*sin(2*pi*y*m/T);
      Rint[:,m+1]   = quadrature.integrate(0,T,numpts,fcos,dim);
      Rint[:,M+m+1] = quadrature.integrate(0,T,numpts,fsin,dim);
    end 
    Rint = Rint/T;
    return Rint;
  end 

  #############################################################
  #  Returns Jacobian of getMHBResidualFC with respect to 
  #  A and T (if varyT == true)
  #############################################################
  function getDMHBResidualFC(A,M,T,mhbpar)
    numpts   = mhbpar.numpts;
    odepar   = mhbpar.odepar;
    dim      = odepar.dim;
    varyT    = mhbpar.varyT;
    if varyT == true
      Jdim = dim*(2*M+1) + 1;
    else
      Jdim = dim*(2*M+1);
    end
    DR       = zeros(dim*(2*M+1),Jdim);
    epsieinv = 1E6;
    epsie    = 1E-6;
    for j = 1 : dim*(2*M+1)
      ej      = zeros(dim*(2*M+1));
      ej[j]   = epsie;
      ej      = reshape(ej,dim,2*M+1);
      DR[:,j] = epsieinv*reshape((getMHBResidualFC(A+ej,M,T,mhbpar)-getMHBResidualFC(A,M,T,mhbpar) ),dim*(2*M+1));
    end
    if varyT == true
      DR[:,end] = epsieinv*reshape((getMHBResidualFC(A,M,T+epsie,mhbpar)-getMHBResidualFC(A,M,T,mhbpar) ),dim*(2*M+1));
    end
    return DR;
  end 

  ######################################################################
  #### Newton iteration for solving getMHBResidualFC = 0 for A or (A,T)
  ######################################################################
  function solveResidualNewton(A,M,T,mhbpar);
    odepar  = mhbpar.odepar;
    dim     = odepar.dim;
    numpts  = mhbpar.numpts;
    tol     = mhbpar.tol;
    maxiter = mhbpar.maxiter;
    varyT    = mhbpar.varyT;
    R    = getMHBResidualFC(A,M,T,mhbpar);
    err  = norm(R,Inf);
    iter = 0;
    while err > tol && iter < maxiter
      DR   = getDMHBResidualFC(A,M,T,mhbpar);
      inc  = -DR \ reshape(R,dim*(2*M+1));
      if varyT == true
        A  = A + reshape(inc[1:end-1],dim,2*M+1);
        T  = T + inc[end];
      else
        A  = A + reshape(inc,dim,2*M+1);
      end
      R    = reshape(getMHBResidualFC(A,M,T,mhbpar),dim*(2*M+1));
      iter = iter + 1;
      err  = norm(R,Inf);
    end 
    if varyT == true
      return (reshape(A,dim,2*M+1),real(T),err,iter);
    else
      return (reshape(A,dim,2*M+1),err,iter);
    end
  end

  struct mhbpar
    odepar ; # An ode par;
    tol    ; # tolerance for the solver
    maxiter; # maximum iterations in the solver
    numpts ; # no. points in the inverse fourier transorm integral
    varyT  ; # true if solving for period T, false otherwise
  end

end
