% This function perform the bisection procedure used in SR-WLS of [14]. 
% And the details of this function is based on reference [4]
% [4] A. Beck, P. Stocia and J. Li, “Exact and approximate
% solutions of source localization problems,” IEEE Trans.
% Signal Process., vol. 56, no. 5, pp. 1770-1778, May 2008.
% [14] S. Tomic, M. Beko, and R. Dinis, "3-D Target Localization
% in Wireless Sensor Networks Using RSS and
% AoA Measurements," in IEEE Transactions on Vehicular
% Technology, vol. 66, no. 4, pp. 3197-3210, April 2017.
% ---------------------------------------------------------
% This code is written bby Kaifeng Gao and Jiang Zhu.
% If you have any problem, please feel free to contact
% jiangzhu16@zju.edu.cn
function y= solve_SRWLS(A,b,D,f,epsilon)
% Input: A,b,D,f: matrix defined in Eq.(17) and Eq.(18) of [4].
%      epsilon: Control algorithm fineness
y_hat = @(lambda) ((A'*A+lambda.*D)^-1)*(A'*b-lambda.*f);
phi = @(lambda) y_hat(lambda)'*D*y_hat(lambda)...
    +2*f'*y_hat(lambda);
% 
[~,ee]=eig(D,A'*A);
low=-(1/ee(end)+1); % 这里要从一个负数开始
% low = -1/(max(diag(ee))+epsilon);
% low=-1;
up=1e3;
lam=low;
while(abs(phi(lam))>epsilon)
    ttt=phi(lam);
    if(ttt>0)
        low=lam;
    else
        up=lam;
    end
    lam=(low+up)/2;
end
y=y_hat(lam);