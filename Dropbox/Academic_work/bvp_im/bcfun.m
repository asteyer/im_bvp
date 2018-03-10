function y = bcfun(x0,xend,d,prob,par)
  y = zeros(d,1);
  if prob == 1 
   % x0(1) = 0 ; xend(1) = 2; 
   y(1) = [1 0]*x0;
   
   y(2) = [1 0]*xend-2;
  end
end

