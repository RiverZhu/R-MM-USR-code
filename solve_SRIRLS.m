function [theta]=solve_SRIRLS(W,A,b,epsilon)
%epsilon=1e-3;
D=[1,   0,  0;
    0,  1,  0;
    0,  0,  0];
f=[0;
   0;
   -0.5];% size(f)==(3,1);
y_hat = @(lambda) ((A'*W*A+lambda.*D)^-1)*(A'*W*b-lambda.*f);
phi = @(lambda) y_hat(lambda)'*D*y_hat(lambda)...
    +2*f'*y_hat(lambda);
% 
[~,ee]=eig(D,A'*W*A);
% low=-1/(ee(end)+epsilon);
low = -1/(max(diag(ee))+epsilon);
% low = max(-diag(A'*W*A));
up=1e3;
lam=low;
while(abs(phi(lam))>epsilon)
    if(phi(lam)>0)
        low=lam;
    else
        up=lam;
    end
    lam=(low+up)/2;
end
% lambda = -6:1e-3:10;
% LL = length(lambda);
% phi_p = zeros(1,LL);
% for i=1:LL
%     phi_p(i) = phi(lambda(i));
% end
% plot(lambda,phi_p)

theta=y_hat(lam);
%phi_lambda=phi(lam);
        
        
        
        
        