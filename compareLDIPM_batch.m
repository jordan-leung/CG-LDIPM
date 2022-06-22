clc
clear all
close all
addpath('/Users/jordan/Google Drive/SCHOOL/PhD/MATLAB/OptimizationFunctions')
pathVar = pathdef;
addpath(pathVar);

% Load data
saveFlag = 0;
caseFlag = 2;
switch caseFlag
    case 1
        load('QPData');
        H = H_QP;
        c = f_QP;
        A = A_QP;
        b = b_QP;
    case 2
        load('QPData2');
        H = H_QP;
        c = f_QP;
        A = A_QP;
        b = b_QP;
    case 3
        N = 10;
        H = diag(linspace(0.01,1,N));
        %         H = eye(N);
        load('minmax')
        xmax = xmax*1000;
        xmin = xmin*1000;
        [A,b] = minmaxMatrices(xmin,xmax);        
end
invH = inv(H);


OPTIONS = optimoptions('quadprog');
OPTIONS = optimoptions(OPTIONS, 'OptimalityTolerance', 1e-10, 'ConstraintTolerance', 1e-10);
x_QP = quadprog(H,c,A,b,[],[],[],[],[],OPTIONS);


% Run LDIPM with normal settings
mu_f_vec = logspace(-12,-2); 
mu_0 = 1e8;
maxIter = 1000;
maxCGIter = 10000;
CGPreCondFlag = 0;
printFlag = 0;
v0 = zeros(size(A,1),1);

numSample = length(mu_f_vec);
execTime_reg = zeros(numSample,1);
execTime_cg = zeros(numSample,1);
execTime_cgdiag = zeros(numSample,1);
NumAverage = 25;
for iOuter = 1:numSample
    execTimeReg_sum = 0;
    execTimeCG_sum = 0;
    execTimeDiag_sum = 0;
    for j = 1:NumAverage
    % Set mu
    mu_f = mu_f_vec(iOuter);
    
    % Reg
    [~,~,~,~,~,~,~,~,~,execTimeReg_i] = logInteriorPoint(H,c,A,b,[],[],mu_f,mu_0,v0,maxIter,printFlag);
    execTimeReg_sum = execTimeReg_sum + execTimeReg_i;
    
    % CGLDIPM
    [~,~,~,~,~,~,~,~,~,execTimeCG_i] = logInteriorPoint_conjgrad_direct(H,c,A,b,mu_f,mu_0,v0,maxIter,maxCGIter,0,printFlag);
    execTimeCG_sum = execTimeCG_sum + execTimeCG_i;
    
    % CGLDIPM with diag
    [~,~,~,~,~,~,~,~,~,execTimeDiag_i] = logInteriorPoint_conjgrad_direct(H,c,A,b,mu_f,mu_0,v0,maxIter,maxCGIter,1,printFlag);
    execTimeDiag_sum = execTimeDiag_sum + execTimeDiag_i;
    end
    execTime_reg(iOuter) = execTimeReg_sum/NumAverage;
    execTime_cg(iOuter) = execTimeCG_sum/NumAverage;
    execTime_cgdiag(iOuter) = execTimeDiag_sum/NumAverage;
end

% dataname = ['.\Data\muPlot_Case',num2str(caseFlag)];
% save(dataname,'mu_f_vec','execTime_reg','execTime_cg')

%% Plotting
close all

figure
h1 = loglog(1./mu_f_vec,execTime_reg,'linewidth',2);
hold on; grid on; box on
h2 = loglog(1./mu_f_vec,execTime_cg,'linewidth',2);
h3 = loglog(1./mu_f_vec,execTime_cgdiag,'linewidth',2);
figSize = [0 0 0.2 0.2];
set(gcf,'units','normalized','position',figSize)
xlabel('$1/\mu_f$','interpreter','latex','fontsize',15)
ylabel('Exec Time (s)','interpreter','latex','fontsize',15)
legend([h1 h2 h3],'Direct','CG','CG Diag','interpreter','latex','fontsize',12,'location','best')
% filename = [pwd '/Figures/muPlot_Case',num2str(caseFlag)];
% saveas(gcf,filename,'epsc')

% figure
% ind = find(mu_f_vec <
% h1 = semilogx(1./mu_f_vec,execTime_reg,'linewidth',2);
% hold on; grid on; box on
% h2 = semilogx(1./mu_f_vec,execTime_cg,'linewidth',2);
% figSize = [0 0 0.2 0.2];
% set(gcf,'units','normalized','position',figSize)
% xlabel('$1/\mu_f$','interpreter','latex','fontsize',15)
% ylabel('Exec Time','interpreter','latex','fontsize',15)
% legend([h1 h2],'Reg','Diag','interpreter','latex','fontsize',12,'location','best')
% ylim([0 Inf])
