function [J,F] = midptbvpsetup(X,T,H,d,prob,par,M)
   J = zeros(d*(M),d*(M)); F=zeros(d*M,1);
    for i= 1:M-1
      % form Jacobian
      Ai =getA((X(i+1,:)+X(i,:))/2,(T(i)+T(i+1))/2,d,prob,par);
      Si = -eye(d)-0.5*H(i)*Ai; Ri = eye(d)-0.5*H(i)*Ai;    
      J((i-1)*d+1:i*d,(i-1)*d+1:i*d) = Si;
      J((i-1)*d+1:i*d,(i)*d+1:(i+1)*d) = Ri;
      % form rhs
      F((i-1)*d+1:i*d) = X(i+1,:)'-X(i,:)'-H(i)*frhs((X(i+1,:)+X(i,:))/2,(T(i)+T(i+1))/2,d,prob,par);
    end
    [B0, BT] = dbcfun(X(1,:)',X(end,:)',d,prob,par);
    J((M-1)*d+1:M*d,1:d) = B0; J((M-1)*d+1:M*d,(M-1)*d+1:M*d) = BT;
    g = bcfun(X(1,:)',X(end,:)',d,prob,par);
    F((M-1)*d+1:end)=g;
end

