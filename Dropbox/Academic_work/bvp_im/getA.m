function A = getA(x,t,d,prob,par )
  A=zeros(d,d);
  if prob == 1
    % y'' + y = 0 , y(0)=0, y(pi/2) = 2, y(x) = 2 sin(x)
    % x = y' implies that x' = y'' = -y
    A = [0 1 ; -1 0];
  end
end

