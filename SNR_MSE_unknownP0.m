% This code is for the following paper:
% Kaifeng Gao, Jiang Zhu and Zhiwei Xu, A Majorization-Minimization based 
% Target Localization Problem from Range Measurements, published 
% in IEEE Communications Letters, 2020
% This code is written bby Kaifeng Gao and Jiang Zhu.
% If you have any problem, please feel free to contact
% jiangzhu16@zju.edu.cn
% ---------------------------------------------------------------
% This simulation compares R-MM-USR with SR-WLS, MSE versus SNR (dB) (Fig.3)
% SR-WLS is performed based on reference [6]
% ---------------------------------------------------------
clear;clc;close all
%% Initialize
x_1 = -2;
x_2 = 3;
x = [x_1;x_2]; % true position of the target, we choose dimension d=2
L0 = 10; % (dB) 
P0 =10^(-L0/10); % Transmit power of the target
gamma = 2;  % path-loss exponent ¦Ã

ai1 = [0;-3;-2;-4;-4];  % Coordinates of sensors.
ai2 = [5;1;2;-1;4];  
aa = [ai1,ai2];
aa = aa';         % size(aa)==(2,N)
m = length(ai1);
r = sqrt(sum((aa-repmat(x,1,m)).^2));
r = r';             % size(r)==(m,1), r=[r_1,r_2,...,r_m]^T
N_obj=200;          % Number of iterations
N_MC=500;           % Number of iterations
SNR=25:2:40;        % Signal to noise ratio (SNR, dB)
L_SNR=length(SNR);
sigma2_n=sum(r.^2)./(m*10.^(SNR/10)); % sigma2_n represents ¦Ò^2
%% Comparison R-MM-USR with SR-WLS for each SNR
% Initialize these matrix to store mean square error (MSE) of x and P0 for
% both approaches.
mse_USR_x=zeros(L_SNR,1);
mse_USR_P0=zeros(L_SNR,1);
mse_SRWLS_x = zeros(L_SNR,1);
mse_SRWLS_P0 = zeros(L_SNR,1);
% Initialize these matrix to store CRB of both x and P0
CRB_x=zeros(L_SNR,1);
CRB_P0=zeros(L_SNR,1);
for nn=1:L_SNR              % loop for SNR
    % calculate FIM
    J=zeros(3,3);
    for i=1:m
        yi=[(x_1-aa(1,i))/r(i);(x_2-aa(2,i))/r(i)];
        J = J + (1/sigma2_n(nn))*...
            [ ( (10*gamma/(r(i)*log(10)))^2 ) *(yi*yi'), (-10/(P0*log(10)))*(10*gamma/(r(i)*log(10))) * yi;
            (-10/(P0*log(10)))*(10*gamma/(r(i)*log(10))) * yi', (10/(P0*log(10)))^2 ];          
    end
    diag_inv_J=diag(J^-1);
    CRB_x(nn) = sum(diag_inv_J(1:2));   % calculate CRB of x
    CRB_P0(nn) = diag_inv_J(end);       % CRB of P0
    
    for mc=1:N_MC    % loop for Monte Carlo trials
        % initialize received signal
        noise=sqrt(sigma2_n(nn))*randn(m,1);
        pL = L0 + 10*gamma*log10(r) + noise;    % path loss model in Eq.(1a) of [14]
        s =  10.^(pL/(10*gamma));               % the corresponding s
                                                % s = [s_1,...,s_m]^T.

        % USR initialize point
        bb=sum(aa.^2,1)';  % bb and AA correspond to the definition in Eq.(24)
        AA=[2*aa',s.^2,repmat(-1,m,1)];
        tilde_x_USR=((AA'*AA)^-1)*AA'*bb;
        x_USR = tilde_x_USR(1:2);
        P0_USR = tilde_x_USR(3);
        
        % iteration for MM algorithm, USR_ini
        x_t=x_USR;
        P0_t = P0_USR;
        for t=1:N_obj
            % fix P0£¬optimize x
            r_t=sqrt(sum((aa-repmat(x_t,1,m)).^2)); % aa=[a_1,a_2,...,a_m] size(r_i)==(d,m)==(2,5)
            r_t = r_t';% size(r_t)==(N,1)
            x_t=[sum( ai1+(s*P0_t^(1/gamma)).*(x_t(1)-ai1)./r_t )/m;
                sum(  ai2+(s*P0_t^(1/gamma)).*(x_t(2)-ai2)./r_t )/m];
            
            % fix x£¬optimize P0
            P0_t = (sum(s.*r_t)/sum(s.^2))^gamma;
        end
        % Calculate mse of x and P0
        mse_USR_x(nn)=mse_USR_x(nn)+sum((x_t-x).^2);
        mse_USR_P0(nn)=mse_USR_P0(nn)+(P0_t-P0).^2;
        
        % Perform SR-WLS according to Section III-B of [14]
        beta = 10.^(-pL/(10*gamma));  % ¦Â=10^(-pL/(10¦Ã));
        % w = 1-pL/sum(pL);
        W_tilde = diag(sqrt(1-pL/sum(pL))); % W = I_3 k-product diag(w);
        A_tilde = [-2*repmat(beta.^2,1,2).*aa',beta.^2,repmat([0,-1],m,1)]; % r0=1m
        A_tilde = W_tilde*A_tilde;
        b_tilde = -(beta.^2).*sum(aa.^2,1)';
        b_tilde = W_tilde*b_tilde;
        small_L_tilde = [0,0,-1/2,0,-1/2]';
        D_tilde = diag([1,1,0,1,0]);
        y_tilde=solve_SRWLS(A_tilde,b_tilde,D_tilde,small_L_tilde,1e-5);
        
        r_hat = sqrt(sum((aa-repmat(y_tilde(1:2),1,m)).^2));
        r_hat = r_hat';
        L0_hat = sum(pL-20*log10(r_hat))/m;
        
        lambda_hat = 10.^((L0_hat - pL)/(10*gamma));
        A1 = [-2*repmat(lambda_hat.^2,1,2).*aa',lambda_hat.^2];
            
        b1 = 1-(lambda_hat.^2).*sum(aa.^2,1)';   
        dij = 10.^((pL-L0_hat)/(10*gamma));
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
mse_USR_x=mse_USR_x/N_MC;
mse_USR_P0=mse_USR_P0/N_MC;
mse_SRWLS_x = mse_SRWLS_x/N_MC;
mse_SRWLS_P0 = mse_SRWLS_P0/N_MC;
%% plot
close all
alw = 0.75;    % AxesLineWidth
fsz = 14;      % Fontsize
lw = 1.5;      % LineWidth
msz = 12;       % MarkerSize

figure;
subplot(2,1,1);  
hold on;box on
set(gca, 'FontSize', fsz, 'LineWidth', alw);
plot(SNR,mse_USR_x,...
    '-ro','LineWidth',lw, 'MarkerSize', msz);
plot(SNR,mse_SRWLS_x,...
    '--b*','LineWidth',lw, 'MarkerSize', msz);
plot(SNR,CRB_x,...
    '-.k+','LineWidth',lw, 'MarkerSize', msz);
legend('R-MM-USR','SR-WLS [6]','CRB (27)');

ylabel('MSE of $\bf{x}$','Interpreter','Latex','Fontsize',fsz)
set(gca,'Yscale','log');
set(gca,'ytick',[1e-4,1e-3,1e-2]);
axis([25,40,5*1e-5,1e-2])

subplot(2,1,2);hold on;box on
set(gca, 'FontSize', fsz, 'LineWidth', alw);
plot(SNR,mse_USR_P0,...
    '-ro','LineWidth',lw, 'MarkerSize', msz);
plot(SNR,mse_SRWLS_P0,...
    '--b*','LineWidth',lw, 'MarkerSize', msz);
plot(SNR,CRB_P0,...
    '-.k+','LineWidth',lw, 'MarkerSize', msz);
xlabel('SNR (dB)');ylabel('\rm{MSE of} $P_0$','Interpreter','Latex','Fontsize',fsz)
set(gca,'Yscale','log');
%% save
% save('mse_snr_uP\mse_USR_x.mat','mse_USR_x');
% save('mse_snr_uP\mse_SRWLS_x.mat','mse_SRWLS_x');
% save('mse_snr_uP\CRB_x.mat','CRB_x');
% 
% save('mse_snr_uP\mse_USR_P0.mat','mse_USR_P0');
% save('mse_snr_uP\mse_SRWLS_P0.mat','mse_SRWLS_P0');
% save('mse_snr_uP\CRB_P0.mat','CRB_P0');
