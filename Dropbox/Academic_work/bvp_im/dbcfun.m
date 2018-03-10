function [dg0,dgend] = dbcfun(x0,xend,d,prob,par)
   if prob == 1
     % x0(1) = 0 ; xend(1) = 2; 
     dg0 = [1 0; 0 0 ]; dgend = [0 0 ; 1 0];
   end
end

