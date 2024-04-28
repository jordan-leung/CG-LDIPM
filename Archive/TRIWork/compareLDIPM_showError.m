clc
clear all
close all

% Load data
saveFlag = 0;
caseFlag = 5;
switch caseFlag
    case 1 % Active constraints MPC example
        load('QPData'); 
        H = H_QP;
        c = f_QP;
        A = A_QP;
        b = b_QP;
    case 2 % Inactive constraints MPC example
        load('QPData2');
        H = H_QP;
        c = f_QP;
        A = A_QP;
        b = b_QP;
    case 3 % Case 1 but with only the input constraints
        load('QPData');
        H = H_QP;
        c = f_QP;
        xmin = -ones(size(H_QP,1),1);
        xmax = ones(size(H_QP,1),1);
        [A,b] = minmaxMatrices(xmin,xmax);
    case 4 % Scalable QP example
        N = 10;
        H = diag(linspace(0.01,1,N));
        %         H = eye(N);
        load('minmax')
        [A,b] = minmaxMatrices(xmin,xmax);
    case 5 % Much more constraints than variables
        
%         % Static variables
%         N = 100;
%         condTarget = 1e1;
%         m = N;
%         A = zeros(m,N);
%         for i = 1:m
%             A(i,1+(i-1)*(N/m):i*(N/m)) = ones(1,N/m);
%         end
%         b_low = -1;
%         b_high = 1;
%         c_low = -10;
%         c_high = 10;
%         
%         % Changing Variables
%         hh= 2*rand(N,N)-1 + 2*rand(N,N)-1;
%         hh = hh*hh';              % symmetric with random entries beween -2 and 2
%         [u, s, v] = svd(hh);
%         s = diag(s);           % s is vector
%         s = s(1)*( 1-((condTarget-1)/condTarget)*(s(1)-s)/(s(1)-s(end))) ;
%         s = diag(s);           % back to matrix
%         H = u*s*v';
%         H = 1/2*(H' + H);
%         c =  c_low + (c_high-c_low).*rand(N,1);
%         b = b_low + (b_high - b_low)*rand(m,1);
        
        %                 save('Case5_Data','H','c','A','b');
        load('Case5_Data');
end
invH = inv(H);


% Run LDIPM with normal settings
mu_f = 1e-10;
% mu_f = 2e-2;
mu_0 = 1e8;
maxIter = 50000;
maxCGIter = 100000;
CGTol = 1e-10;
CGPreCondFlag = 0;
printFlag = 1;
v0 = zeros(size(A,1),1);
if caseFlag < 5
    vThresh = -4;
    vNumThresh = length(v0)/4;
else
    vThresh = 0;
    vNumThresh = -1; % so always on
end

% Given parameters (delta,epsilon), deteremine the shortstep parameters (N,k)
m = size(A,1);
delta = .45; % delta in (0,1/2)
epsilon = delta; % epsilon in (0,invq(delta)), invq(x) and x intersect at like 0.9 so this is always valid
zeta = acosh(delta/2 + 1) - epsilon;
if epsilon == delta
    N_ls = 1;
else
    N_ls = ceil(1 + log(epsilon)/(log(delta)*log(2)));
end
k_ls = exp(2*acosh(zeta^2/(2*m) +1 ));
beta = 1/2;
params.delta = delta;
params.gamma = delta/2 + sqrt(delta^2/4 + delta);
params.epsilon = epsilon;
params.zeta = zeta;
params.N = N_ls;
params.k = k_ls;
params.beta = beta;
const.params = params;

% Longstep CG
[x_ls,lambda_ls,s_ls,v_ls,mu_ls,startFlag_ls,numIter_ls,muStar_ls,exitFlag_ls,execTime_ls,CGIters_ls,CGres_ls,warmStartRes_ls,CGerror_ls] = logInteriorPoint_conjgrad(H,c,A,b,mu_f,mu_0,v0,maxIter,maxCGIter,CGTol,1,1,vThresh,vNumThresh);

% % Shortstep_optimal - no WS
% fprintf('------ Running with CGLDIPM-Shortstep (Optimal)  ------ \n')
% [x,v,mu,execTime,numIter,CGIters,CGres,CGerror,dHist,dDiffHist,dInitHist] = logInteriorPoint_conjgrad_shortStep_opt(H,c,A,b,mu_f,mu_0,v0,maxIter,maxCGIter,CGTol,1,params,vThresh,vNumThresh,1);

% PID Mu
fprintf('------ Running with CGLDIPM-Shortstep (D Invariant)  ------ \n') %  W,c,Aineq,bineq,mu_f,mu_0,v0,maxIter,maxCGIter,CGTol,preCondFlag,beta,vThresh,vNumThresh,wsFlag)
dInvar = 1;
[x1,v1,mu1,execTime1,numIter1,CGIters1,CGres1,CGerror1,dHist1,dDiffHist1,dInitHist1,stepHist] = logInteriorPoint_dInvariant(H,c,A,b,mu_f,mu_0,v0,maxIter,maxCGIter,CGTol,1,dInvar,vThresh,vNumThresh,1);

%% Plotting
close all


% --------------- Plots for the shortstep method ---------------

% figure
% h1 = plot(CGIters,'linewidth',2);
% grid on; box on; hold on;
% h2 = plot(CGIters1,'linewidth',2);
% xlabel('LDIPM Iteration','interpreter','latex','fontsize',15)
% ylabel('CG Iterations','interpreter','latex','fontsize',15)
% legend('Shortstep','New','interpreter','latex','fontsize',12,'location','best')
% ylim([1 Inf])
% 
% figure
% semilogy(dHist,'linewidth',2);
% grid on; box on; hold on;
% semilogy(dHist1,'linewidth',2);
% xlabel('LDIPM Iteration','interpreter','latex','fontsize',15)
% ylabel('$\|d\|_2$','interpreter','latex','fontsize',15)
% legend('Shortstep','New','interpreter','latex','fontsize',12,'location','best')
% 
% figure
% semilogy(mu,'linewidth',2);
% grid on; box on; hold on;
% semilogy(mu1,'linewidth',2);
% xlabel('LDIPM Iteration','interpreter','latex','fontsize',15)
% ylabel('$\mu$','interpreter','latex','fontsize',15)
% legend('Shortstep','New','interpreter','latex','fontsize',12,'location','best')
% 
% figure
% semilogy(mu,'linewidth',2);
% grid on; box on; hold on;
% semilogy(mu1,'linewidth',2);
% xlabel('LDIPM Iteration','interpreter','latex','fontsize',15)
% ylabel('$\mu$','interpreter','latex','fontsize',15)
% legend('Shortstep','New','interpreter','latex','fontsize',12,'location','best')
% xlim([1 600])

sum(sum(CGIters_ls))
sum(CGIters1)

% figure
% set(gcf,'units','normalized','position',[0.6 0.5 0.1 0.2])
% 
% figure
% set(gcf,'units','normalized','position',[0.6 0.5 0.1 0.2])
% 
% figure
% set(gcf,'units','normalized','position',[0.6 0.5 0.1 0.2])
% 
% figure
% set(gcf,'units','normalized','position',[0.6 0.5 0.1 0.2])
% figure
% set(gcf,'units','normalized','position',[0.6 0.5 0.1 0.2])
% figure
% set(gcf,'units','normalized','position',[0.6 0.5 0.1 0.2])



figure
h2 = plot(CGIters1,'k','linewidth',2);
grid on; box on; hold on;
xlabel('LDIPM Iteration','interpreter','latex','fontsize',15)
ylabel('CG Iterations','interpreter','latex','fontsize',15)
ylim([1 Inf])
xlim([1 Inf])
set(gcf,'units','normalized','position',[0.6 0.5 0.1 0.2])

figure
subplot(2,1,1)
semilogy(dHist1(:,1),'k','linewidth',2);
grid on; box on; hold on;
ylabel('$\|d\|_2$','interpreter','latex','fontsize',15)
subplot(2,1,2);
semilogy(dHist1(:,2),'k','linewidth',2);
box on; grid on;
xlabel('LDIPM Iteration','interpreter','latex','fontsize',15)
ylabel('$\|d\|_\infty$','interpreter','latex','fontsize',15)
xlim([1 Inf])
set(gcf,'units','normalized','position',[0.6 0.5 0.1 0.2])

figure
semilogy(mu1,'k','linewidth',2);
grid on; box on; hold on;
xlabel('LDIPM Iteration','interpreter','latex','fontsize',15)
ylabel('$\mu$','interpreter','latex','fontsize',15)
set(gcf,'units','normalized','position',[0.6 0.5 0.1 0.2])


figure
plot(1./(stepHist(:,1)+stepHist(:,2)),'k','linewidth',2);
grid on; box on; hold on;
xlabel('LDIPM Iteration','interpreter','latex','fontsize',15)
ylabel('$c$','interpreter','latex','fontsize',15)
plot(1./(stepHist(:,1)),'--','linewidth',2);
plot(1./(1+stepHist(:,2)),'-.','linewidth',2);
legend('$\frac{1}{a+b}$','$\frac{1}{a}$','$\frac{1}{1+b}$','interpreter','latex','fontsize',15,'location','best')
xlim([1 Inf])
set(gcf,'units','normalized','position',[0.6 0.5 0.1 0.2])

















