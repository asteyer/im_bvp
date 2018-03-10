function y = frhs(x,t,d,par,prob)
  y = zeros(d,1);
  if prob == 1
      % y'' + y = 0 , y(0)=0, y(pi/2) = 2, y(x) = 2 sin(x)
      % x = y' implies that x' = y'' = -y
      y(1) = x(2);
      y(2) = -x(1);
  end

end


