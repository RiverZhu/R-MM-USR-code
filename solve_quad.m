function theta = solve_quad(N,bb,ee,cc,epsilon)
x_hat = @(lambda) (-1/(N+lambda))*(bb+lambda*ee); 
phi = @(lambda) x_hat(lambda)'*x_hat(lambda)...
    +2*ee'*x_hat(lambda)+cc;
low=-N;
up=1e20;
lam=low;
while(abs(phi(lam))>epsilon)
    if(phi(lam)>0)
        low=lam;
    else
        up=lam;
    end
    lam=(low+up)/2;
end
theta=x_hat(lam);

