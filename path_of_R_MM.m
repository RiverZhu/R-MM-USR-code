% This code is for the following paper:
% Kaifeng Gao, Jiang Zhu and Zhiwei Xu, A Majorization-Minimization based 
% Target Localization Problem from Range Measurements, published 
% in IEEE Communications Letters, 2020
% This code is written bby Kaifeng Gao and Jiang Zhu.
% If you have any problem, please feel free to contact
% jiangzhu16@zju.edu.cn
% ---------------------------------------------------------------
% This simulation compares different initial points of R-MM (Fig.1 of the paper).
% ---------------------------------------------------------
clear;clc;close all;
%% Initialize
x_1 = -2;
x_2 = 3;
x = [x_1;x_2];% true position of the target, we choose dimension d=2

ai1=[6;0;5;1;3];
ai2=[4;-10;-3;-4;-3];
m=length(ai1); % number of sensors
aa = [ai1,ai2];
aa = aa';         % size(r_i)==(2,N)
r = sqrt(sum((aa-repmat(x,1,m)).^2));
r = r';             % size(r)==(m,1), r=[r_1,r_2,...,r_m]^T
SNR=40;
sigma2_n=sum(r.^2)./(m*10.^(SNR/10));
%%
% initialize received signal
noise=sqrt(sigma2_n)*randn(m,1);
s=r+noise;
% USR-LS initialize point
A=[-2*aa',ones(m,1)]; %size(A)==(N,3)£¬for USR initialize point
N_obj=100;
b=s.^2-sum(aa.^2)';
y_star=((A'*A)^-1)*A'*b;
x_USR=y_star(1:2); % USR initial point
x_ini1 = [6;-2]; % The initial point of case 1
x_ini2 = [4;-2]; % The initial point of case 2
%%
% Initialize these matrix to store the points of the iteration path
path_USR = zeros(2,N_obj);
path_ini1 = zeros(2,N_obj);
path_ini2 = zeros(2,N_obj);
% Iteration of USR
x_t=x_USR;
for t=1:N_obj
    path_USR(:,t)=x_t;
    r_t=sqrt(sum((aa-repmat(x_t,1,m)).^2));
    r_t = r_t';% size(d_USR)==(N,1)
    x_t=[sum(ai1+s.*(x_t(1)-ai1)./r_t)/m;
        sum(ai2+s.*(x_t(2)-ai2)./r_t)/m];
end
% Iteration of case 1
x_t=x_ini1;
for t=1:N_obj
    path_ini1(:,t)=x_t;
    r_t=sqrt(sum((aa-repmat(x_t,1,m)).^2));
    r_t = r_t';% size(d_USR)==(N,1)
    x_t=[sum(ai1+s.*(x_t(1)-ai1)./r_t)/m;
        sum(ai2+s.*(x_t(2)-ai2)./r_t)/m];
end

% Iteration of constraint case 1
x_t=x_ini2;
for t=1:N_obj
    path_ini2(:,t)=x_t;
    r_t=sqrt(sum((aa-repmat(x_t,1,m)).^2));
    r_t = r_t';% size(d_USR)==(N,1)
    x_t=[sum(ai1+s.*(x_t(1)-ai1)./r_t)/m;
        sum(ai2+s.*(x_t(2)-ai2)./r_t)/m];
end
% Pick a few points in the iteration path, otherwise it is too dense.
index_USR=[1,2:4:20];
USR_p_x=path_USR(1,index_USR);
USR_p_y=path_USR(2,index_USR);

index1 = [1:4,5:3:15,18:5:28];
ini1_p_x1=path_ini1(1,index1); % x_1 for ploting of  case 1
ini1_p_x2=path_ini1(2,index1); % x_2 for ploting of  case 1

index2 = [1:8,10:4:24];
ini2_p_x1=path_ini2(1,index2);
ini2_p_x2=path_ini2(2,index2);
%% Calculate the contour of objective function for ploting
obj_fun = @(r_i,s,theta) sum((s-sqrt(sum((r_i-repmat(theta,1,m)).^2))').^2);
x1_scan=linspace(-5,12,1000);
x2_scan=linspace(-6,6,1000);
Lx=length(x1_scan);
Ly=length(x2_scan);
fun_p=zeros(Ly,Lx);
for ii=1:Lx
    for jj=1:Ly
        fun_p(jj,ii)=obj_fun(aa,s,[x1_scan(ii);x2_scan(jj)]);
    end
end
%%
close all
alw = 0.75;    % AxesLineWidth
fsz = 14;      % Fontsize
lw = 1.4;      % LineWidth
msz = 8;      % MarkerSize

figure; hold on;box on
set(gca, 'FontSize', fsz, 'LineWidth', alw);

% plot path_USR, path_1 and path_2 
plot(ini1_p_x1,ini1_p_x2,'-r+','LineWidth',lw,'MarkerSize', msz);
plot(ini2_p_x1,ini2_p_x2,':k*','LineWidth',lw,'MarkerSize', msz);
plot(USR_p_x,USR_p_y,':m+','LineWidth',lw,'MarkerSize', msz);

% plot global optimum and local optimum
plot(x_1,x_2,'ro','LineWidth',lw,'MarkerSize', msz);
plot(path_ini1(1,end),path_ini1(2,end),'ks','LineWidth',lw,'MarkerSize', msz); 
% plot the contour of objective function
contour(x1_scan,x2_scan,fun_p,20);
xlabel('$x_1$','Interpreter','Latex','Fontsize',fsz)
ylabel('$x_2$','Interpreter','Latex','Fontsize',fsz)
legend('R-MM-\bf{x}_1^{(0)}','R-MM-\bf{x}_2^{(0)}','R-MM-USR','Global Optimum','Local Optimum');
axis equal
