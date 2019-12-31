% This code is for the following paper:
% Kaifeng Gao, Jiang Zhu and Zhiwei Xu, A Majorization-Minimization based 
% Target Localization Problem from Range Measurements, published 
% in IEEE Communications Letters, 2020
% This code is written bby Kaifeng Gao and Jiang Zhu.
% If you have any problem, please feel free to contact
% jiangzhu16@zju.edu.cn
% ---------------------------------------------------------------
% This simulation compares R-MM-USR with SR-WLS when changing ¦Ã (Fig.4).
% SR-WLS is performed based on reference [6]
% ---------------------------------------------------------
clear;clc;close all
%% Initialize
x_1 = -2; % true position of the target
x_2 = 3;  % we choose the dimension d=2
x = [x_1;x_2]; 

L0 = 10; % (dB) Path-loss value at a short reference distance r_0, 
         % which in the paper we set it as 1 m.
P0 =10^(-L0/10); % Transmit power of the target, L0 = -10 log10(P0)

ai1 = [3;1;-1;-4;-4];  % Coordinates of sensors.
ai2 = [5;6;2;-1;4];   
m = length(ai1);       % Number of sensors.
aa = [ai1,ai2];        % sort the coordinates into a matrix.
aa = aa';              % size(aa)==(2,N)
r = sqrt(sum((aa-repmat(x,1,m)).^2));
r = r';                % size(r)==(m,1), r=[r_1,r_2,...,r_m]^T
N_obj=400;             % Number of iterations
N_MC=200;              % Number of Monte Carlo trials
SNR=40;                % Signal to noise ratio (SNR, dB)
sigma2_n=sum(r.^2)./(m*10.^(SNR/10)); % sigma2_n represents ¦Ò^2
gamma=1.8:0.1:2.2;
L_gamma=length(gamma);
%% Comparison R-MM-USR with SR-WLS for each gamma
% Initialize these matrix to store mean square error (MSE) of x and P0 for
% both approaches.
mse_USR_x = zeros(L_gamma,1);
mse_USR_P0 = zeros(L_gamma,1);
mse_SRWLS_x = zeros(L_gamma,1);
mse_SRWLS_P0 = zeros(L_gamma,1);
for nn=1:L_gamma  % loop for gamma
    for mc=1:N_MC % loop for Monte Carlo trials
        
        % initialize received signal
        noise=sqrt(sigma2_n)*randn(m,1);
        pL = L0 + 10*gamma(nn)*log10(r) + noise; % path loss model in Eq.(1a) of [14]
        s =  10.^(pL/(10*gamma(nn)));            % the corresponding s
                                                 % s = [s_1,...,s_m]^T.

        % USR initialize point
        bb=sum(aa.^2,1)';               % bb and AA correspond to the 
        AA=[2*aa',s.^2,repmat(-1,m,1)]; % definition in Eq.(24)
        tilde_x_USR=((AA'*AA)^-1)*AA'*bb;
        x_USR = tilde_x_USR(1:2);
        P0_USR = tilde_x_USR(3);

        % iteration for MM algorithm 
        x_t=x_USR;
        P0_t = P0_USR;
        for t=1:N_obj
            % fix P0£¬optimize x
            r_t=sqrt(sum((aa-repmat(x_t,1,m)).^2)); 
            r_t = r_t';% size(r_t)==(N,1)
            % update x 
            x_t=[sum( ai1+(s*P0_t^(1/gamma(nn))).*(x_t(1)-ai1)./r_t )/m;
                sum(  ai2+(s*P0_t^(1/gamma(nn))).*(x_t(2)-ai2)./r_t )/m];
            
            % fix x£¬optimize P0
            P0_t = (sum(s.*r_t)/sum(s.^2))^gamma(nn); % update P0 as Eq.(23)
        end
        % Calculate mse of x and P0
        mse_USR_x(nn)=mse_USR_x(nn)+sum((x_t-x).^2);
        mse_USR_P0(nn)=mse_USR_P0(nn)+(P0_t-P0).^2;
 
        % Perform SR-WLS according to Section III-B of [14].      
        % 1) Solve Eq.(15) of [14] to obtain the initial estimate of x
        beta = 10.^(-pL/(10*gamma(nn)));  % ¦Â=10^(-pL/(10¦Ã));
        % w_tilde = 1-pL/sum(pL); % W_tilde = I_3 Kronecker-product diag(w_tilde);
        W_tilde = diag(sqrt(1-pL/sum(pL))); 
        A_tilde = [-2*repmat(beta.^2,1,2).*aa',beta.^2,repmat([0,-1],m,1)]; % r0=1m
        A_tilde = W_tilde*A_tilde;
        b_tilde = -(beta.^2).*sum(aa.^2,1)';
        b_tilde = W_tilde*b_tilde;
        small_L_tilde = [0,0,-1/2,0,-1/2]';
        D_tilde = diag([1,1,0,1,0]);
        y_tilde=solve_SRWLS(A_tilde,b_tilde,D_tilde,small_L_tilde,1e-5);
        
        % 2) Use x_hat (r_hat) to compute the ML estimate of L0
        r_hat = sqrt(sum((aa-repmat(y_tilde(1:2),1,m)).^2));
        r_hat = r_hat';
        L0_hat = sum(pL-20*log10(r_hat))/m;
        
        % 3) ExploitL0 to calculate lambda_hat and  use this estimated
        % value to solve Eq.(12) of [14].
        lambda_hat = 10.^((L0_hat - pL)/(10*gamma(nn)));
        A1 = [-2*repmat(lambda_hat.^2,1,2).*aa',lambda_hat.^2];    
        b1 = 1-(lambda_hat.^2).*sum(aa.^2,1)';   
        dij = 10.^((pL-L0_hat)/(10*gamma(nn)));
        W1 = diag(1-dij/sum(dij));
        A1 = W1*A1;
        b1 = W1*b1;
        D1 = diag([1,1,0]); small_L1=[0,0,-1/2]';
        y1= solve_SRWLS(A1,b1,D1,small_L1,1e-5);
        % Calculate mse of x and P0
        x_SRWLS=y1(1:2);
        mse_SRWLS_x(nn) = mse_SRWLS_x(nn)+sum((x_SRWLS-x).^2);
        P0_hat = 10^(-L0_hat/10);
        mse_SRWLS_P0(nn) = mse_SRWLS_P0(nn) + (P0_hat - P0).^2; 
    end
end
% Calculate the average of the results of Monte Carlo
mse_USR_x = mse_USR_x/N_MC;
mse_USR_P0 = mse_USR_P0/N_MC;
mse_SRWLS_x = mse_SRWLS_x/N_MC;
mse_SRWLS_P0 = mse_SRWLS_P0/N_MC;
%% plot
close all
alw = 0.75;    % AxesLineWidth
fsz = 14;      % Fontsize
lw = 1.5;      % LineWidth
msz = 12;       % MarkerSize

figure;
subplot(2,1,1);hold on;box on
set(gca, 'FontSize', fsz, 'LineWidth', alw);
plot(gamma,mse_USR_x,...
    '-ro','LineWidth',lw, 'MarkerSize', msz);
plot(gamma,mse_SRWLS_x,...
    '--b*','LineWidth',lw, 'MarkerSize', msz);
legend('R-MM-USR','SR-WLS [6]');
ylabel('MSE of $\bf{x}$','Interpreter','Latex','Fontsize',fsz)
set(gca,'Yscale','log');
axis([gamma(1),gamma(end),1e-4,1e-1]);

subplot(2,1,2);hold on;box on
set(gca, 'FontSize', fsz, 'LineWidth', alw);
plot(gamma,mse_USR_P0,...
    '-ro','LineWidth',lw, 'MarkerSize', msz);
plot(gamma,mse_SRWLS_P0,...
    '--b*','LineWidth',lw, 'MarkerSize', msz);
xlabel('$\gamma$','Interpreter','Latex','Fontsize',fsz');
ylabel('\rm{MSE of} $P_0$','Interpreter','Latex','Fontsize',fsz)
set(gca,'Yscale','log');
axis([gamma(1),gamma(end),5*1e-8,1e-3]);
set(gca,'ytick',[1e-7,1e-5,1e-3]);
%% save
% save('mse_snr_uP\mse_USR_x.mat','mse_USR_x');
% save('mse_snr_uP\mse_SRWLS_x.mat','mse_SRWLS_x');
% save('mse_snr_uP\CRB_x.mat','CRB_x');
% 
% save('mse_snr_uP\mse_USR_P0.mat','mse_USR_P0');
% save('mse_snr_uP\mse_SRWLS_P0.mat','mse_SRWLS_P0');
% save('mse_snr_uP\CRB_P0.mat','CRB_P0');
