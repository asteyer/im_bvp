module quadrature  ### this module implements a fastish version of high order gauss quadrature in 1D

  #################################################3
  ## This integrates a function f:R-> R^dim
  ## over the interval [a,b] with a quadrature rule 
  ## defined by quad_rule using an invterval divided 
  ## into numquadpts subintervals
  ##
  function integrate(a,b,numquadpts,f,dim)
    dx = (b-a)/numquadpts;
    Rint = zeros(dim);
    x1 = a; x2 = a + dx;
    for j = 1:numquadpts
      Rint = Rint + quad_rule(x1,x2,f,dim)*dx;
      x1 = x2;
      x2 = x1 + dx;
    end
    return Rint;
  end

  ################ 
  # quad_rule integrates f:R->R^dim on the interval [a,b] 
  # with a higher order Gauss quadrature rule
  ################
  function quad_rule(a,b,f,dim)
    fbar = xbar->f(0.5*(b-a)*xbar + 0.5*(a+b));
    Rint = (1/10)*fbar(-1)[:] + (49/90)*fbar(-sqrt(3/7))[:] + (32/45)*fbar(0)[:] + (49/90)*fbar(sqrt(3/7))[:] + (1/10)*fbar(1)[:]
    return Rint;
  end 

end
