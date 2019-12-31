function [theta]=solve_SRLS(A,b,epsilon)
%epsilon=1e-3;
D=[1,   0,  0;
    0,  1,  0;
    0,  0,  0];
f=[0;
   0;
   -0.5];% size(f)==(3,1);
y_hat = @(lambda) ((A'*A+lambda.*D)^-1)*(A'*b-lambda.*f);
phi = @(lambda) y_hat(lambda)'*D*y_hat(lambda)...
    +2*f'*y_hat(lambda);
% 
[~,ee]=eig(D,A'*A);
% low=-1/ee(end)+1e-6;
low = -1/(max(diag(ee))+epsilon);
up=1e10;
lam=low;
while(abs(phi(lam))>epsilon)
    if(phi(lam)>0)
        low=lam;
    else
        up=lam;
    end
    lam=(low+up)/2;
end
theta=y_hat(lam);
theta=theta(1:2);
%phi_lambda=phi(lam);
        
        
        
        
        