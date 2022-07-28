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
delta = .01; % delta in (0,1/2)
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

% Shortstep_optimal - no WS
fprintf('------ Running with CGLDIPM-Shortstep (Optimal)  ------ \n')
[x,v,mu,execTime,numIter,CGIters,CGres,CGerror,dHist,dDiffHist,dInitHist] = logInteriorPoint_conjgrad_shortStep_opt(H,c,A,b,mu_f,mu_0,v0,maxIter,maxCGIter,CGTol,1,params,vThresh,vNumThresh,1);

 % Regular
fprintf('------ Running with CGLDIPM-Shortstep  ------ \n')
[x2,v2,mu2,execTime2,numIter2,CGIters2,CGres2,CGerror2,dHist2] = logInteriorPoint_conjgrad_shortStep(H,c,A,b,mu_f,mu_0,v0,maxIter,maxCGIter,CGTol,1,params,vThresh,vNumThresh,1);

% Regular criterion but with explicit error
fprintf('------ Running with CGLDIPM-Shortstep (Calc Error) ------ \n')
[x3,v3,mu3,execTime3,numIter3,CGIters3,CGres3,CGerror3,dHist3] = logInteriorPoint_conjgrad_shortStep_calcError(H,c,A,b,mu_f,mu_0,v0,maxIter,maxCGIter,CGTol,1,params,vThresh,vNumThresh);
% 
% % Regular criterion but with explicit divergence
% fprintf('------ Running with CGLDIPM-Shortstep (Calc Error) ------ \n')
% [x4,v4,mu4,execTime4,numIter4,CGIters4,CGres4,CGerror4,dHist4] = logInteriorPoint_conjgrad_shortStep_explicit_onlyDiv(H,c,A,b,mu_f,mu_0,v0,maxIter,maxCGIter,CGTol,1,params,vThresh,vNumThresh);
% 
%  % Regular criterion but with explicit divergence and error
% fprintf('------ Running with CGLDIPM Calc Divergence  ------ \n')
% [x5,v5,mu5,execTime5,numIter5,CGIters5,CGres5,CGerror5,dHist5,hHist5] = logInteriorPoint_conjgrad_shortStep_explicit(H,c,A,b,mu_f,mu_0,v0,maxIter,maxCGIter,CGTol,1,params,vThresh,vNumThresh);
% 
% % % Other criterion with explicit divergence (compare to onlyDiv)
% % fprintf('------ Running with CGLDIPM-Shortstep (Calc Error) ------ \n')
% % [x6,v6,mu6,execTime6,numIter6,CGIters6,CGres6,CGerror6,dHist6] = logInteriorPoint_conjgrad_shortStep_explicit2(H,c,A,b,mu_f,mu_0,v0,maxIter,maxCGIter,CGTol,1,params,vThresh,vNumThresh);

% % % Using additional M-norm based bound on e
% % fprintf('------ Running with CGLDIPM-Shortstep (Other bound) ------ \n')
% [x7,v7,mu7,execTime7,numIter7,CGIters7,CGres7,CGerror7,dHist7] = logInteriorPoint_conjgrad_shortStep_secondErrorBound(H,c,A,b,mu_f,mu_0,v0,maxIter,maxCGIter,CGTol,1,params,vThresh,vNumThresh);



% fixed = sum(CGIters)
% early = sum(CGIters2)
% errorbased = sum(CGIters3)
% divbased = sum(CGIters4)

%% Plotting
close all


% --------------- Plots for the shortstep method ---------------

figure
set(gcf,'units','normalized','position',[0.1300 0.1192 0.1 0.2])
h1 = plot(CGIters,'linewidth',2);
grid on; box on; hold on;
h2 = plot(CGIters2,'linewidth',2);
h3 = plot(CGIters3,'linewidth',2);
h4 = plot(CGIters4,'linewidth',2);
h5 = plot(CGIters5,'linewidth',2);
xlabel('LDIPM Iteration','interpreter','latex','fontsize',15)
ylabel('CG Iterations','interpreter','latex','fontsize',15)
legend('Fixed truncation','Regular Bound','Bound w/ Exact Error','Bound w/ Exact Div','Bound w/ Exact Error and Div','interpreter','latex','fontsize',12,'location','best')
ylim([1 Inf])

figure
set(gcf,'units','normalized','position',[0.1300 0.1192 0.1 0.2])
semilogy(CGres,'linewidth',2)
grid on; box on; hold on;
semilogy(CGres2,'linewidth',2);
semilogy(CGres3,'linewidth',2);
semilogy(CGres4,'linewidth',2);
semilogy(CGres5,'linewidth',2);
xlabel('LDIPM Iteration','interpreter','latex','fontsize',15)
ylabel('CG Residual','interpreter','latex','fontsize',15)
legend('Fixed truncation','Regular Bound','Bound w/ Exact Error','Bound w/ Exact Div','Bound w/ Exact Error and Div','interpreter','latex','fontsize',12,'location','best')
ylim([1e-12 Inf])


figure
set(gcf,'units','normalized','position',[0.1300 0.1192 0.1 0.2])
semilogy(CGerror,'linewidth',2)
grid on; box on; hold on;
semilogy(CGerror2,'linewidth',2);
semilogy(CGerror3,'linewidth',2);
semilogy(CGerror4,'linewidth',2);
semilogy(CGerror5,'linewidth',2);
xlabel('LDIPM Iteration','interpreter','latex','fontsize',15)
ylabel('CG Error','interpreter','latex','fontsize',15)
legend('Fixed truncation','Regular Bound','Bound w/ Exact Error','Bound w/ Exact Div','Bound w/ Exact Error and Div','interpreter','latex','fontsize',12,'location','best')
ylim([1e-12 Inf])


% --------------- Compare the other bound ---------------

% figure
% set(gcf,'units','normalized','position',[0.1300 0.1192 0.1 0.2])
% h1 = plot(CGIters,'linewidth',2);
% grid on; box on; hold on;
% plot(CGIters4,'linewidth',2);
% plot(CGIters6,'linewidth',2);
% xlabel('LDIPM Iteration','interpreter','latex','fontsize',15)
% ylabel('CG Iterations','interpreter','latex','fontsize',15)
% legend('Fixed truncation','Bound 1 w/ Exact Div','Bound 2 w/ Exact Div','interpreter','latex','fontsize',12,'location','best')
% ylim([1 Inf])
% 
% figure
% set(gcf,'units','normalized','position',[0.1300 0.1192 0.1 0.2])
% semilogy(CGres,'linewidth',2)
% grid on; box on; hold on;
% semilogy(CGres4,'linewidth',2);
% semilogy(CGres6,'linewidth',2);
% xlabel('LDIPM Iteration','interpreter','latex','fontsize',15)
% ylabel('CG Residual','interpreter','latex','fontsize',15)
% legend('Fixed truncation','Bound 1 w/ Exact Div','Bound 2 w/ Exact Div','interpreter','latex','fontsize',12,'location','best')
% ylim([1e-12 Inf])


% % --------------- Compare the bounds on e ---------------
% figure
% set(gcf,'units','normalized','position',[0.1300 0.1192 0.1 0.2])
% h1 = plot(CGIters,'linewidth',2);
% grid on; box on; hold on;
% h2 = plot(CGIters2,'linewidth',2);
% h5 = plot(CGIters7,'linewidth',2);
% xlabel('LDIPM Iteration','interpreter','latex','fontsize',15)
% ylabel('CG Iterations','interpreter','latex','fontsize',15)
% legend('Fixed truncation','Residual Bound','Residual and M-norm Bound','interpreter','latex','fontsize',12,'location','best')
% ylim([1 Inf])
% 
% figure
% set(gcf,'units','normalized','position',[0.1300 0.1192 0.1 0.2])
% semilogy(CGres,'linewidth',2)
% grid on; box on; hold on;
% semilogy(CGres2,'linewidth',2);
% semilogy(CGres7,'linewidth',2);
% xlabel('LDIPM Iteration','interpreter','latex','fontsize',15)
% ylabel('CG Residual','interpreter','latex','fontsize',15)
% legend('Fixed truncation','Residual Bound','Residual and M-norm Bound','interpreter','latex','fontsize',12,'location','best')
% ylim([1e-12 Inf])
% 
% 
% figure
% set(gcf,'units','normalized','position',[0.1300 0.1192 0.1 0.2])
% semilogy(CGerror,'linewidth',2)
% grid on; box on; hold on;
% semilogy(CGerror2,'linewidth',2);
% semilogy(CGerror7,'linewidth',2);
% xlabel('LDIPM Iteration','interpreter','latex','fontsize',15)
% ylabel('CG Error','interpreter','latex','fontsize',15)
% legend('Fixed truncation','Residual Bound','Residual and M-norm Bound','interpreter','latex','fontsize',12,'location','best')
% ylim([1e-12 Inf])
% 
% 
% figure
% set(gcf,'units','normalized','position',[0.1300 0.1192 0.1 0.2])
% semilogy(dHist,'linewidth',2)
% grid on; box on; hold on;
% semilogy(dHist2,'linewidth',2);
% semilogy(dHist7,'linewidth',2);
% plot([1 length(dHist)],[params.gamma params.gamma])
% xlabel('LDIPM Iteration','interpreter','latex','fontsize',15)
% ylabel('$\|d\|$','interpreter','latex','fontsize',15)
% legend('Fixed truncation','Residual Bound','Residual and M-norm Bound','$\|d\|$ upper-bound','interpreter','latex','fontsize',12,'location','best')
% xlim([1 Inf])

% -------------- For d diff --------------------------
figure
indVec = find(~isinf(dDiffHist));
set(gcf,'units','normalized','position',[0.1300 0.1192 0.1 0.2])
semilogy(dDiffHist(indVec),'linewidth',2);
grid on; box on; hold on;
xlabel('LDIPM Iteration','interpreter','latex','fontsize',15)
ylabel('$\| d_{i+1} - d_i \|$','interpreter','latex','fontsize',15)



