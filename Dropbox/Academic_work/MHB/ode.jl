module ode

  function frhs(state,time,par)
    problem   = par.problem;
    dim       = par.dim;
    modelpars = par.modelpars;
    statedot  = zeros(dim);
    if problem == 1 #### Linear, harmonic oscillator 
      w = modelpars[1]
      statedot = [0 w ; -w 0]*state;
    elseif problem == 2  #### Example from page 162 of Hirsch, Smale, Devaney
      statedot[1] = 0.5*state[1] - state[2] - 0.5 * ( state[1]^3 + state[1]*state[2]^2);
      statedot[2] = state[1] + 0.5*state[2] - 0.5 * ( state[2]^3 + state[2]*state[1]^2);
    end
    return statedot;
  end
  
  function Dfrhs(state,tim,par)
    problem   = par.problem;
    dim       = par.dim;
    modelpars = par.modelpars;
    if problem == 1
      J = [0 1; -1 0];
    end
    return J;
  end

  struct odepar
    dim; # Int dimension of state
    problem; # Int specifying which problem we are solving
    modelpars; # vector of model pars
  end

end
