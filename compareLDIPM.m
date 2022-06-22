clc
clear all
close all
addpath('/Users/jordan/Google Drive/SCHOOL/PhD/MATLAB/OptimizationFunctions')
pathVar = pathdef;
addpath(pathVar);

% Load data
saveFlag = 0;
caseFlag = 1;
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
    case 3 % Scalable QP example
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
x_QP = quadprog(H,c,A,b,[],[],[],[],[],OPTIONS)


% Run LDIPM with normal settings
mu_f = 1e-10;
mu_0 = 1e8;
maxIter = 200;
maxCGIter = 10000;
CGPreCondFlag = 0;
printFlag = 1;
v0 = zeros(size(A,1),1);

fprintf('------ Running with regular scheme ------ \n')
[xReg,lambdaReg,sReg,vReg,~,~,numIterReg,~,~,execTimeReg] = logInteriorPoint(H,c,A,b,[],[],mu_f,mu_0,v0,maxIter,printFlag);

% % Diagonal noPrecond
% fprintf('------ Running with CGLDIPM (no precond) ------ \n')
% [x,~,~,vStar,muStar,~,numIter,~,~,execTime,CGIters,CGres,CGpchist,wsRes] = logInteriorPoint_conjgrad(H,c,A,b,mu_f,mu_0,v0,maxIter,maxCGIter,0,printFlag);
% 
% % Diagonal precond
% fprintf('------ Running with CGLDIPM (diag precond) ------ \n')
% [x2,~,~,vStar2,muStar2,~,numIter2,~,~,execTime2,CGIters2,CGres2,CGpchist2,wsRes2] = logInteriorPoint_conjgrad(H,c,A,b,mu_f,mu_0,v0,maxIter,maxCGIter,1,printFlag);


% % Diagonal noPrecond
% fprintf('------ Running with CGLDIPM (no precond) ------ \n')
% [x,~,~,vStar,muStar,~,numIter,~,~,execTime,CGIters,CGres,CGpchist,wsRes] = logInteriorPoint_conjgrad_pproj(H,c,A,b,mu_f,mu_0,v0,maxIter,maxCGIter,0,printFlag);
% 
% % Diagonal precond
% fprintf('------ Running with CGLDIPM (diag precond) ------ \n')
% [x2,~,~,vStar2,muStar2,~,numIter2,~,~,execTime2,CGIters2,CGres2,CGpchist2,wsRes2] = logInteriorPoint_conjgrad_pproj(H,c,A,b,mu_f,mu_0,v0,maxIter,maxCGIter,1,printFlag);

% Diagonal noPrecond
fprintf('------ Running with CGLDIPM (no precond) ------ \n')
[x,~,~,vStar,muStar,~,numIter,~,~,execTime,CGIters,CGres,CGpchist,wsRes] = logInteriorPoint_conjgrad_pproj_doubleWS(H,c,A,b,mu_f,mu_0,v0,maxIter,maxCGIter,0,printFlag);

% Diagonal precond
fprintf('------ Running with CGLDIPM (diag precond) ------ \n')
[x2,~,~,vStar2,muStar2,~,numIter2,~,~,execTime2,CGIters2,CGres2,CGpchist2,wsRes2] = logInteriorPoint_conjgrad_pproj_doubleWS(H,c,A,b,mu_f,mu_0,v0,maxIter,maxCGIter,1,printFlag);



%% Plotting
close all

% Warm-start comparison plot
figure
markSize = 15;
hold on; grid on; box on
h1 = plot(CGIters(:,1),'linewidth',2);
h2 = plot(CGIters(:,2),'Color',h1.Color,'linestyle','-.','linewidth',2);
h3 = plot(CGIters2(:,1),'linewidth',2);
h4 = plot(CGIters2(:,2),'Color',h3.Color,'linestyle','-.','linewidth',2);
plot(CGIters(:,1),'.','markersize',markSize,'color',h1.Color);
plot(CGIters(:,2),'.','markersize',markSize,'color',h2.Color);
plot(CGIters2(:,1),'.','markersize',markSize,'color',h3.Color);
plot(CGIters2(:,2),'.','markersize',markSize,'color',h4.Color);
figSize = [0 0 0.2 0.2];
set(gcf,'units','normalized','position',figSize)
xlabel('LDIPM Iteration','interpreter','latex','fontsize',15)
ylabel('\# of CG Iterations','interpreter','latex','fontsize',15)
% legend([h1 h2 h3 h4],'No cond, 1st','No cond, 2nd','Diag 1st','Diag 2nd','interpreter','latex','fontsize',12,'location','best')
% ylim([0 18])


figure
markSize = 15;
h5 = semilogy(wsRes,'linewidth',2,'color',h1.Color);
hold on; grid on; box on
h6 = semilogy(wsRes2,'linewidth',2,'color',h3.Color);
semilogy(wsRes,'.','markersize',markSize,'color',h5.Color);
semilogy(wsRes2,'.','markersize',markSize,'color',h6.Color);
figSize = [0 0 0.2 0.2];
set(gcf,'units','normalized','position',figSize)
xlabel('LDIPM Iteration','interpreter','latex','fontsize',15)
ylabel('Warm-start residual','interpreter','latex','fontsize',15)
% legend([h5 h6],'No cond','Diag','interpreter','latex','fontsize',12,'location','best')


% figure
% hold on; grid on; box on
% h1 = plot(CGIters,'.','Markersize',15);
% h2 = plot(CGIters2,'.','Markersize',15);
% % h3 = plot(CGIters3,'.','Markersize',15);
% h4 = plot(CGIters4,'.','Markersize',15);
% plot(CGIters,'Color',h1.Color)
% plot(CGIters2,'Color',h2.Color)
% % plot(CGIters3,'Color',h3.Color)
% plot(CGIters4,'Color',h4.Color)
% figSize = [0 0 0.2 0.2];
% set(gcf,'units','normalized','position',figSize)
% xlabel('LDIPM Iteration','interpreter','latex','fontsize',15)
% ylabel('\# of CG Iterations','interpreter','latex','fontsize',15)
% legend([h1 h2 h4],'None','Diag','Block','interpreter','latex','fontsize',12,'location','best')
% % legend([h1 h2],'None','Diag','interpreter','latex','fontsize',12,'location','best')
% ylim([0 Inf])


