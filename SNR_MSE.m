% This simulation compare the MSE of R-MM-USR,USR-LS, SDR, SR-LS, vs. SNR
% This code is for the following paper:
% Kaifeng Gao, Jiang Zhu and Zhiwei Xu, A Majorization-Minimization based 
% Target Localization Problem from Range Measurements, published 
% in IEEE Communications Letters, 2020
% This code is written bby Kaifeng Gao and Jiang Zhu.
% If you have any problem, please feel free to contact
% jiangzhu16@zju.edu.cn
% ---------------------------------------------------------------
% This simulation compares the MSE of R-MM-USR, USR-LS, SDR, 
% SR-LS, SR-IRLS, versus SNR (dB) (Fig.2 of the paper)
% SR-LS is performed based on reference [5]
% SR-IRLS is performed based on reference [7]
% ---------------------------------------------------------
clear;clc;close all;
%% Initialize
x_1 = -2;
x_2 = 3;
x = [x_1;x_2];
ai1=[6;0;5;1;3];
ai2=[4;-10;-3;-4;-3];
m=length(ai1);
aa = [ai1,ai2];
aa = aa';               % size(aa)==(2,m)
r = sqrt(sum((aa-repmat(x,1,m)).^2));
r = r';                 % size(r)==(m,1),r=[r_1,r_2,...,r_m]^T
A=[-2*aa',ones(m,1)];   % size(A)==(m,3)£¬for USR-LS initialize point
N_obj=200;              % Number of iterations
N_MC=100;               % Number of iterations
SNR=10:10:50;           % Signal to noise ratio (SNR, dB)
sigma2_n=sum(r.^2)./(m*10.^(SNR/10)); % sigma2_n represents ¦Ò^2
L_SNR=length(SNR);
%%
% Initialize these matrix to store mean square error (MSE)
mse_MM_USR=zeros(L_SNR,1);
mse_USRLS=zeros(L_SNR,1);
mse_SRLS=zeros(L_SNR,1);
mse_SDR=zeros(L_SNR,1);
mse_SR_IRLS = zeros(L_SNR,1);
CRB=zeros(L_SNR,1);
% warning off  % avoid warning when using CVX tool box
for nn=1:L_SNR % loop for SNR
    % calculate CRB
    J=zeros(2,2);
    for i=1:m
        yi=[(x_1-aa(1,i))/r(i);(x_2-aa(2,i))/r(i)];
        J = J +(1/sigma2_n(nn))*(yi*yi');
    end
    CRB(nn)=trace(J^-1);
    
    for mc=1:N_MC   % % loop for Monte Carlo trials
        % initialize received signal
        noise=sqrt(sigma2_n(nn))*randn(m,1);
        s=r+noise;
        
        % USR-LS 
        b=s.^2-sum(aa.^2,1)';
        y_star=((A'*A)^-1)*A'*b;
        x_USR=y_star(1:2);
        mse_USRLS(nn)=mse_USRLS(nn)+sum((x_USR-x).^2);
        
        % iteration for R-MM-USR algorithm
        x_t=x_USR;
        for t=1:N_obj
            r_t=sqrt(sum((aa-repmat(x_t,1,m)).^2));
            r_t = r_t';% size(r_t)==(N,1)
            x_t=[sum(ai1+s.*(x_t(1)-ai1)./r_t)/m;
                sum(ai2+s.*(x_t(2)-ai2)./r_t)/m];
        end
        mse_MM_USR(nn)=mse_MM_USR(nn)+sum((x_t-x).^2);
       
        % The inaccurate solutions to the R-LS problem produced by SDR
        % we load the matrix in the following, which has been calculated.
        % you can re-calculate SDR by uncommenting the following two lines.
%         [theta_cvx,~]=cxv_obj(aa,s,m);
%         mse_SDR(nn)=mse_SDR(nn)+sum((theta_cvx-theta).^2);

        % the exact solution of SR-LS
        x_SRLS=solve_SRLS(A,b,1e-8);
        mse_SRLS(nn)=mse_SRLS(nn)+sum((x_SRLS-x).^2);
        % iteration for SR-IRLS (iteratively reweighted least square)
        % we load the matrix in the following, which has been calculated.
        % you can re-calculate SR-IRLS by uncommenting the following lines (line 84~108).
%         W_t = eye(m);
%         y_last = ones(3,1);    
%         y_t = zeros(3,1);
%         count = 0;
%         while sum(abs(y_t - y_last)) > 1e-8 
%             y_last = y_t;
%             y_t = solve_SRIRLS(W_t,A,b,1e-8);
%             for i = 1:m
%                 e_it = [-2*aa(:,i)',1]*y_t - (s(i)^2 - aa(:,i)'*aa(:,i));
%                 epsilon_ = 1.34*sqrt(3*sigma2_n(nn));
%                 w_it = 1/(e_it^2+ epsilon_^2);
%                 W_t(i,i) = w_it;
%             end
%             count = count + 1;
%             if count > 500
%                 fprintf('break');
%                 break;
%             end
%         end 
%         x_SR_IRLS = y_t(1:2);
%         mse_SR_IRLS(nn) = mse_SR_IRLS(nn) + sum((x_SR_IRLS - x).^2);
%         fprintf('x_USR = [%.4f,%.4f]^T   ',x_USR);
%         fprintf('x_SR_IRLS = [%.4f,%.4f]^T count = %d \n',x_SR_IRLS,count);
%         fprintf('nn=%d,L_SNR = %d,mc = %d, N_MC = %d, all: %d / %d \n',...
%             nn,L_SNR,mc,N_MC,mc+(nn-1)*N_MC,N_MC*L_SNR);
    end
end
% warning on  % avoid warning when using CVX tool box
mse_USRLS=mse_USRLS/N_MC;
mse_MM_USR=mse_MM_USR/N_MC;
mse_SRLS=mse_SRLS/N_MC;
mse_SR_IRLS = mse_SR_IRLS/N_MC;
mse_SDR_p=mse_SDR/N_MC;
load('mse_SDR_p.mat');   % if re-calculate SDR, uncomment this line.
load('mse_SR_IRLS.mat'); % if re-calculate SR-IRLS, uncomment this line.
%% plot
close all
alw = 0.75;    % AxesLineWidth
fsz = 14;      % Fontsize
lw = 1.5;      % LineWidth
msz = 12;       % MarkerSize

figure;hold on;box on
set(gca, 'FontSize', fsz, 'LineWidth', alw);
plot(SNR,mse_MM_USR,...
    '-ro','LineWidth',lw, 'MarkerSize', msz);
plot(SNR,mse_USRLS,...
    '-r+','LineWidth',lw, 'MarkerSize', msz);
plot(SNR,mse_SRLS,...
    ':b*','LineWidth',lw, 'MarkerSize', msz);
plot(SNR,mse_SR_IRLS,...
    ':ms','LineWidth',lw, 'MarkerSize', msz);
plot(SNR,mse_SDR_p,...
    '--ko','LineWidth',lw, 'MarkerSize', msz);
plot(SNR,CRB,...
    '-.k+','LineWidth',lw, 'MarkerSize', msz);
legend('R-MM-USR','USR-LS [5]','SR-LS [5]','SR-IRLS [7]','SDR [5]','CRB (6)');
xlabel('SNR (dB)');ylabel('MSE');
set(gca,'Yscale','log');
axis([10,50,1e-3,1e3])
set(gca,'ytick',[1e-3,1e-2,1e-1,1,1e1,1e2,1e3]);
%% save 
save('mse_SR_IRLS.mat','mse_SR_IRLS');