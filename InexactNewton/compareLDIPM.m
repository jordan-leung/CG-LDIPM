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
        H = H_QP;CG
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
        load('Case5_Data');
end
invH = inv(H);

% LDIPM settings
mu_f = 1e-8;
mu_0 = 1e8;
maxIter = 50000;
maxCGIter = 100000;
CGTol = 1e-8;
v0 = zeros(size(A,1),1);

% fprintf('------ Running with regular LDIPM  ------ \n')
% [x_reg,~,~,v_reg] = logInteriorPoint_test(H,c,A,b,mu_f,mu_0,v0,maxIter,1);

fprintf('------ Running with CG LDIPM (longstep)  ------ \n')
[x_cg,~,~,v_cg,muVec_cg,~,~,~,~,~,CGIters_cg,CGres_cg] = logInteriorPoint_conjgrad(H,c,A,b,mu_f,mu_0,v0,maxIter,maxCGIter,CGTol,0,1,1,1);

fprintf('------ Running with CGLDIPM-GIN  ------ \n')
[x,lambda,s,v,muVec,startFlag,numIter,muStar,exitFlag,execTime,CGIters,CGres,FVec,feasVec] = logInteriorPoint_conjgrad_INB(H,c,A,b,mu_f,mu_0,v0,maxIter,maxCGIter,CGTol,1,1000);
% 
fprintf('------ Running with CGLDIPM-GIN  ------ \n')
[x2,~,~,v2,muVec2,~,numIter,~,~,~,CGIters2,CGres2,FVec2,feasVec2] = logInteriorPoint_conjgrad_INB(H,c,A,b,mu_f,mu_0,v0,maxIter,maxCGIter,CGTol,1,900);
% 
fprintf('------ Running with CGLDIPM-GIN  ------ \n')
[x3,~,~,v3,muVec3,~,numIter,~,~,~,CGIters3,CGres3,FVec3,feasVec3] = logInteriorPoint_conjgrad_INB(H,c,A,b,mu_f,mu_0,v0,maxIter,maxCGIter,CGTol,1,800);




%% PLotting
close all
set(0,'defaultAxesFontSize',12)
figSize = [0 0 0.2 0.2];
saveFigFlag = 0;

% Mu vs. CG iterations 
figure
totalIter_cg = cumsum((sum(CGIters_cg'))');
totalIter = cumsum(CGIters);
totalIter2 = cumsum(CGIters2);
totalIter3 = cumsum(CGIters3);
set(gcf,'units','normalized','position',figSize)
h1 = semilogy(totalIter_cg,muVec_cg);
hold on
h2 = semilogy(totalIter,muVec);
h3 = semilogy(totalIter2,muVec2);
h4 = semilogy(totalIter3,muVec3);
semilogy(totalIter_cg,muVec_cg,'.','Markersize',15,'Color',h1.Color)
semilogy(totalIter,muVec,'.','Markersize',15,'Color',h2.Color)
semilogy(totalIter2,muVec2,'.','Markersize',15,'Color',h3.Color)
semilogy(totalIter3,muVec3,'.','Markersize',15,'Color',h4.Color)
semilogy(totalIter_cg(end),muVec_cg(end),'s','Markersize',15,'Linewidth',2,'Color',h1.Color)
semilogy(totalIter(end),muVec(end),'s','Markersize',15,'Linewidth',2,'Color',h2.Color)
semilogy(totalIter2(end),muVec2(end),'s','Markersize',15,'Linewidth',2,'Color',h3.Color)
semilogy(totalIter3(end),muVec3(end),'s','Markersize',15,'Linewidth',2,'Color',h4.Color)
grid on; box on;
hold off
legend([h1 h2 h3 h4],'Longstep w/ CG','Inexact Newton, $\epsilon = 10$',...
    'Inexact Newton, $\epsilon = 1$','Inexact Newton, $\epsilon = 0.5$',...
    'interpreter','latex','fontsize',12,'location','northeast')
xlabel('CG Iterations','interpreter','latex','fontsize',15)
ylabel('$\mu$','interpreter','latex','fontsize',15)
if saveFigFlag
    filename = strcat('./Figures/','preliminaryPlot');
    saveas(gcf,filename,'epsc'); 
end

% CG iterations vs. LDIPM iteration
set(0, 'DefaultLineLineWidth', 2);
figure
CGIters_longstep =  CGIters_cg(:,1) + CGIters_cg(:,2);
set(gcf,'units','normalized','position',figSize)
plot(CGIters_longstep,'color',h1.Color);
grid on; box on; hold on
plot(CGIters,'color',h2.Color);
plot(CGIters2,'color',h3.Color);
plot(CGIters3,'color',h4.Color);
plot(length(CGIters_longstep),CGIters_longstep(end),'s','Markersize',15,'Linewidth',2,'Color',h1.Color)
plot(length(CGIters),CGIters(end),'s','Markersize',15,'Linewidth',2,'Color',h2.Color)
plot(length(CGIters2),CGIters2(end),'s','Markersize',15,'Linewidth',2,'Color',h3.Color)
plot(length(CGIters3),CGIters3(end),'s','Markersize',15,'Linewidth',2,'Color',h4.Color)
legend('Longstep w/ CG','Inexact Newton, $\epsilon = 10$',...
    'Inexact Newton, $\epsilon = 1$','Inexact Newton, $\epsilon = 0.5$',...
    'interpreter','latex','fontsize',12,'location','northeast')
xlabel('LDIPM Iterations','interpreter','latex','fontsize',15)
ylabel('CG Iteratios','interpreter','latex','fontsize',15)
if saveFigFlag
    filename = strcat('./Figures/','preliminaryPlot_iters');
    saveas(gcf,filename,'epsc'); 
end

% F bound
figure
set(gcf,'units','normalized','position',figSize)
plot(feasVec,'color',h2.Color);
hold on; box on; grid on;
plot(feasVec2,'color',h3.Color);
plot(feasVec3,'color',h4.Color);
plot(length(feasVec),feasVec(end),'s','Markersize',15,'Linewidth',2,'Color',h2.Color)
plot(length(feasVec2),feasVec2(end),'s','Markersize',15,'Linewidth',2,'Color',h3.Color)
plot(length(feasVec3),feasVec3(end),'s','Markersize',15,'Linewidth',2,'Color',h4.Color)
xlabel('LDIPM Iterations','interpreter','latex','fontsize',15)
ylabel('Primal Feasibility','interpreter','latex','fontsize',15)
if saveFigFlag
    filename = strcat('./Figures/','preliminaryPlot_feas');
    saveas(gcf,filename,'epsc'); 
end


