clc
clear all
close all

% Load data
saveFlag = 0;
caseFlag = 3;
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

% % First, solve for the unconstrained solution and project onto solution
% xStar_unc = -H\c;
% xProj = xStar_unc;
% n = size(H,1);
% m = size(A,1);
% for i = 1:n
%    if xProj(i) >  xmax(i)
%        xProj(i) = xmax(i);
%    elseif xProj(i) < xmin(i)
%        xProj(i) = xmin(i);
%    end
% end
% 
% % Calculate slack
% s_init = A*xProj + b;
% for i = 1:m
%     if s_init(i) < 1e-6
%         s_init(i) = 1e-6;
%     end
% end
% v0 = -log(s_init);

fprintf('------ Running with regular scheme ------ \n')
[xReg,lambdaReg,sReg,vReg,~,~,numIterReg,~,~,execTimeReg] = logInteriorPoint(H,c,A,b,[],[],mu_f,mu_0,v0,maxIter,printFlag);

% Diagonal noPrecond
fprintf('------ Running with CGLDIPM (no precond) ------ \n')
[x,~,~,vStar,muStar,~,numIter,~,~,execTime,CGIters,CGres,wsRes,CGerror] = logInteriorPoint_conjgrad(H,c,A,b,mu_f,mu_0,v0,maxIter,maxCGIter,0,printFlag);

% Diagonal precond
fprintf('------ Running with CGLDIPM (diag precond) ------ \n')
[x2,~,~,vStar2,muStar2,~,numIter2,~,~,execTime2,CGIters2,CGres2,wsRes2,CGerror] = logInteriorPoint_conjgrad(H,c,A,b,mu_f,mu_0,v0,maxIter,maxCGIter,1,printFlag);

if caseFlag == 3
    [x3,iterCount3,xError3,execTime3] = projGradSolver_rt(H,c,zeros(size(H,1)),xmin,xmax,1000,x_QP,1e-6);
end

%% Plotting
close all

% Warm-start comparison plot
figure
markSize = 15;
hold on; grid on; box on
h1 = plot(CGIters(:,1),'linewidth',2);
h2 = plot(CGIters(:,2),'Color',h1.Color,'linestyle','-.','linewidth',2);
h3 = plot(CGIters2(:,1),'linewidth',2);
h4 = plot(CGIters2(:,2),'Color',h3.Color,'linestyle','--','linewidth',2);
plot(CGIters(:,1),'.','markersize',markSize,'color',h1.Color);
plot(CGIters(:,2),'.','markersize',markSize,'color',h2.Color);
plot(CGIters2(:,1),'.','markersize',markSize,'color',h3.Color);
plot(CGIters2(:,2),'.','markersize',markSize,'color',h4.Color);
figSize = [0 0 0.2 0.2];
set(gcf,'units','normalized','position',figSize)
xlabel('LDIPM Iteration','interpreter','latex','fontsize',15)
ylabel('\# of CG Iterations','interpreter','latex','fontsize',15)
% legend([h1 h2 h3 h4],'No cond, 1st','No cond, 2nd','Diag 1st','Diag 2nd','interpreter','latex','fontsize',12,'location','best')
ylim([0 20])

figure
markSize = 15;
hold on; grid on; box on
h5 = plot(CGIters(:,1)+CGIters(:,2),'linewidth',2,'color',h1.Color);
h6 = plot(CGIters2(:,1)+CGIters2(:,2),'linewidth',2,'color',h3.Color);
plot(CGIters(:,1)+CGIters(:,2),'.','markersize',markSize,'color',h1.Color);
plot(CGIters2(:,1)+CGIters2(:,2),'.','markersize',markSize,'color',h3.Color);
figSize = [0 0 0.2 0.2];
set(gcf,'units','normalized','position',figSize)
xlabel('LDIPM Iteration','interpreter','latex','fontsize',15)
ylabel('Total CG Iterations','interpreter','latex','fontsize',15)
% legend([h5 h6],'No cond','Diag','interpreter','latex','fontsize',12,'location','best')
ylim([0 Inf])


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


% % Error plot
% figure
% markSize = 15;
% h1 = semilogy(CGerror(:,1),'linewidth',2);
% hold on; grid on; box on
% h2 = semilogy(CGerror(:,2),'Color',h1.Color,'linestyle','-.','linewidth',2);
% h3 = semilogy(CGerror2(:,1),'linewidth',2);
% h4 = semilogy(CGerror2(:,2),'Color',h3.Color,'linestyle','--','linewidth',2);
% semilogy(CGerror(:,1),'.','markersize',markSize,'color',h1.Color);
% semilogy(CGerror(:,2),'.','markersize',markSize,'color',h2.Color);
% semilogy(CGerror2(:,1),'.','markersize',markSize,'color',h3.Color);
% semilogy(CGerror2(:,2),'.','markersize',markSize,'color',h4.Color);
% semilogy = [0 0 0.2 0.2];
% set(gcf,'units','normalized','position',figSize)
% xlabel('LDIPM Iteration','interpreter','latex','fontsize',15)
% ylabel('CG Error','interpreter','latex','fontsize',15)
% ylim([0 20])


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


