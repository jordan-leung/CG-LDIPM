clc
clear all
close all

% Load data
saveFigFlag = 0;
caseFlag = 6;
switch caseFlag
    case 1 % Active constraints MPC example
        load('QPData'); 
        H = H_QP;
        c = f_QP;
        A = A_QP;
        b = b_QP;
    case 2 % Scalable QP example
        N = 10;
        H = diag(linspace(0.01,1,N));
        %         H = eye(N);
        load('minmax')
        [A,b] = minmaxMatrices(xmin,xmax);
    case 3 % Much more constraints than variables
        load('Case5_Data');
    case 4
        load('randData') % n = 200, m = 500
    case 5
        load('randData2') % n = 50, m = 150
    case 6
        load('MPCData')
        H = H_QP;
        c = f_QP;
        A = A_QP;
        b = b_QP;
end
invH = inv(H);

% LDIPM settings
mu_f = 1e-4;
mu_0 = 1e6;
maxIter = 150;
maxCGIter = 1000;
CGTol = 1e-6;
v0 = zeros(size(A,1),1);

% Set options
opts.mu_f = mu_f;
opts.mu_0 = mu_0;
opts.maxIter = maxIter;
opts.printFlag = 1;
opts.maxCGIter = maxCGIter;
opts.CGTol = CGTol;

% Set shortstep parameters
m  = size(A,1);
gamma = 0.9;
phi = 0.99*sqrt(1-4*(1-gamma)*gamma^(-2));
theta = phi^2*gamma^2/(2 - gamma + sqrt(4*(1-gamma) + phi^2*gamma^2));
epsilon  = 0.5*qInv(theta);
N = ceil(log(epsilon)/log(theta));
zeta = qInv(theta) - epsilon;
kappa = exp(2*qInv(zeta^2/m));
opts.N = N;
opts.kappa = kappa;
opts.gamma = gamma;
opts.zeta = zeta;

% Get centered point
fprintf('------ Acquiring centred point  ------ \n')
opts_init = opts;
opts_init.mu_f = mu_0;
dTol =  1e-4;
[x_c,v_c,d] = logInteriorPoint_getCenteredPoint(H,c,A,b,v0,mu_0,1000,dTol);
if norm(d,'inf') > dTol
    error('Didnt get initial centred point')
end

% Reconfigurate to use a manual kappa
opts.N = 1;
opts.kappa = (1/0.2);

% Get centered point using full CG
% fprintf('------ Running longstep LDIPM  ------ \n')
% [xStar,outputStar] = logInteriorPoint(H,c,A,b,v_c,opts);
% fprintf('------ Running optimal CG-LDIPM  ------ \n')
% [x1,output1] = cgLDIPM(H,c,A,b,v_c,opts);
% fprintf('------ Running inexact CG-LDIPM  ------ \n')
% [x2,output2] = cgLDIPM_inexact(H,c,A,b,v_c,opts);
% 
fprintf('------ Running longstep LDIPM  ------ \n')
[xStar,outputStar] = cgLDIPM_longstep(H,c,A,b,v_c,opts);
fprintf('------ Running exact medium step CG-LDIPM  ------ \n')
[x1,output1] = cgLDIPM_medStep(H,c,A,b,v_c,opts);
fprintf('------ Running inexact medium step CG-LDIPM  ------ \n')
[x2,output2] = cgLDIPM_medStep_inexact(H,c,A,b,v_c,opts);
% fprintf('------ Running exact new step CG-LDIPM  ------ \n')
% [x3,output3] = cgLDIPM_newStep(H,c,A,b,v_c,opts);

norm(x1-xStar)
norm(x2-xStar)
% norm(x3-xStar)

%% Plotting
close all

% Post-process to get plot of mu (y-axis) vs. total cg iterations
CGIters = output1.CGIters; % for base
totalIters = sum(CGIters);
muVecTotal1 = zeros(totalIters,1);
count = 1;
for i = 1:length(CGIters)
    for j = 1:CGIters(i)
        muVecTotal1(count) = output1.muVec(i);
        count = count + 1;
    end
end
CGIters = output2.CGIters; % for inexact
totalIters = sum(CGIters);
muVecTotal2 = zeros(totalIters,1);
count = 1;
for i = 1:length(CGIters)
    for j = 1:CGIters(i)
        muVecTotal2(count) = output2.muVec(i);
        count = count + 1;
    end
end
CGIters = outputStar.CGIters; % for longstep
totalIters = sum(CGIters);
muVecTotal3 = zeros(totalIters,1);
count = 1;
for i = 1:length(CGIters)
    for j = 1:CGIters(i)
        muVecTotal3(count) = outputStar.muVec(i);
        count = count + 1;
    end
end

% Plotting
set(0,'defaultLineLineWidth',1)
set(0,'defaultAxesFontSize',12)
labelsize = 16;
legendsize = 12;
colorMatrix = colororder;
colorMatrix(7,:) = [0 102 0]/255;
subPlotGap = [0.07 0.1];
subPlotH = [0.1 0.05];
subPlotW = [0.1 0.05];
figSize = [0 0 0.3 0.5];

figure
set(gcf,'units','normalized','position',figSize)
subtightplot(4,1,1,subPlotGap,subPlotH,subPlotW)
h1 = semilogy(muVecTotal1,'color',colorMatrix(1,:));
box on; grid on; hold on
semilogy(length(muVecTotal1),muVecTotal1(end),'.' ,'color',colorMatrix(1,:),'markersize',18);
h2 = semilogy(muVecTotal2,'color',colorMatrix(2,:));
semilogy(length(muVecTotal2),muVecTotal2(end),'.' ,'color',colorMatrix(2,:),'markersize',18);
h3 = semilogy(muVecTotal3,'color',colorMatrix(3,:));
semilogy(length(muVecTotal3),muVecTotal3(end),'.' ,'color',colorMatrix(3,:),'markersize',18);
xlabel('Cumulative CG Iterations','interpreter','Latex','FontSize',labelsize)
ylabel('$\mu$','interpreter','Latex','FontSize',labelsize)
legend([h1 h2 h3],'Nominal','Inexact','Longstep','interpreter','Latex','FontSize',legendsize)
yticks([1e-4 1e0 1e4])
xlim([1 Inf])

subtightplot(4,1,2,subPlotGap,subPlotH,subPlotW)
h2 = semilogy([1 length(output2.condVec(:,1))],[CGTol CGTol],'linewidth',1,'color',colorMatrix(1,:));
box on; grid on; hold on
h3 = semilogy(output2.condVec(:,1),'linewidth',1,'color',colorMatrix(2,:));
% h4 = semilogy(output2.condVec(:,2),'linewidth',1,'color',colorMatrix(3,:));
xlabel('LDIPM Iteration','interpreter','Latex','FontSize',labelsize)
ylabel('Residual','interpreter','Latex','FontSize',labelsize)
legend('Nominal Requirement','Inexact Requirement','interpreter','Latex','FontSize',legendsize,'location','southwest')
% legend('Base Requirement','Requirement 1','Requirement 2','interpreter','Latex','FontSize',legendsize,'location','southwest')
ylim([1e-9 Inf])
xlim([1 Inf])
switch caseFlag
    case 3
        yticks([1e-8 1e-4 1e0])
    case 4
        yticks([1e-8 1e-4 1e0])
    case 6
        yticks([1e-8 1e-4 1e0])
end

subtightplot(4,2,5,subPlotGap,subPlotH,subPlotW)
semilogy(output1.muVec,'-','color',colorMatrix(1,:));
box on; grid on; hold on
semilogy(output2.muVec,'color',colorMatrix(2,:));
semilogy(outputStar.muVec,'color',colorMatrix(3,:));
xlabel('LDIPM Iteration','interpreter','Latex','FontSize',labelsize)
ylabel('$\mu$','interpreter','Latex','FontSize',labelsize)
xlim([1 Inf])


subtightplot(4,2,6,subPlotGap,subPlotH,subPlotW)
plot(output1.CGIters,'-','color',colorMatrix(1,:));
box on; grid on; hold on
plot(output2.CGIters,'color',colorMatrix(2,:));
plot(outputStar.CGIters,'color',colorMatrix(3,:));
xlabel('LDIPM Iteration','interpreter','Latex','FontSize',labelsize)
ylabel('CG Iterations','interpreter','Latex','FontSize',labelsize)
xlim([1 Inf])




subtightplot(4,2,7,subPlotGap,subPlotH,subPlotW)
semilogy(output2.sigma,'color',colorMatrix(2,:))
box on; grid on; hold on
xlabel('LDIPM Iteration','interpreter','Latex','FontSize',labelsize)
ylabel('$\sigma_{\mu}(x)$','interpreter','Latex','FontSize',labelsize)
switch caseFlag
    case 3
        yticks([1e-10 1e-5 1e0])
    case 4
        yticks([1e-4 1e-2 1e0 1e2])
    case 6
        yticks([1e-4 1e0 1e4 1e8])
end
xlim([1 Inf])

subtightplot(4,2,8,subPlotGap,subPlotH,subPlotW)
semilogy(output2.minCond,'color',colorMatrix(2,:))
box on; grid on; hold on
xlabel('LDIPM Iteration','interpreter','Latex','FontSize',labelsize)
ylabel('$\delta(\|d\|)$','interpreter','Latex','FontSize',labelsize)
xlim([1 Inf])

% if saveFigFlag == 1
%     filename = strcat('./Plots/shortstep',num2str(caseFlag));
%     saveas(gcf,filename,'epsc')
% end

if saveFigFlag == 1
    filename = strcat('./Plots/comparison',num2str(caseFlag));
    saveas(gcf,filename,'epsc')
end

%% Plotting for TRI pres
close all

%
% Post-process to get plot of mu (y-axis) vs. total cg iterations
CGIters = output1.CGIters; % for base
totalIters = sum(CGIters);
muVecTotal1 = zeros(totalIters,1);
count = 1;
for i = 1:length(CGIters)
    for j = 1:CGIters(i)
        muVecTotal1(count) = output1.muVec(i);
        count = count + 1;
    end
end
CGIters = output2.CGIters; % for inexact
totalIters = sum(CGIters);
muVecTotal2 = zeros(totalIters,1);
count = 1;
for i = 1:length(CGIters)
    for j = 1:CGIters(i)
        muVecTotal2(count) = output2.muVec(i);
        count = count + 1;
    end
end
CGIters = outputStar.CGIters; % for longstep
totalIters = sum(CGIters);
muVecTotal3 = zeros(totalIters,1);
count = 1;
for i = 1:length(CGIters)
    for j = 1:CGIters(i)
        muVecTotal3(count) = outputStar.muVec(i);
        count = count + 1;
    end
end

% Plotting
set(0,'defaultLineLineWidth',1)
set(0,'defaultAxesFontSize',12)
labelsize = 16;
legendsize = 12;
colorMatrix = colororder;
colorMatrix(7,:) = [0 102 0]/255;
subPlotGap = [0.09 0.1];
subPlotH = [0.1 0.05];
subPlotW = [0.1 0.05];
figSize = [0 0 0.3 0.4];

figure
set(gcf,'units','normalized','position',figSize)
subtightplot(3,1,1,subPlotGap,subPlotH,subPlotW)
h1 = semilogy(muVecTotal1,'color',colorMatrix(1,:));
box on; grid on; hold on
semilogy(length(muVecTotal1),muVecTotal1(end),'.' ,'color',colorMatrix(1,:),'markersize',18);
h2 = semilogy(muVecTotal2,'color',colorMatrix(2,:));
semilogy(length(muVecTotal2),muVecTotal2(end),'.' ,'color',colorMatrix(2,:),'markersize',18);
h3 = semilogy(muVecTotal3,'color',colorMatrix(3,:));
semilogy(length(muVecTotal3),muVecTotal3(end),'.' ,'color',colorMatrix(3,:),'markersize',18);
xlabel('Cumulative CG Iterations','interpreter','Latex','FontSize',labelsize)
ylabel('$\mu$','interpreter','Latex','FontSize',labelsize)
legend([h1 h2 h3],'Nominal','Inexact','Longstep','interpreter','Latex','FontSize',legendsize)
yticks([1e-4 1e0 1e4])
xlim([1 Inf])

subtightplot(3,1,2,subPlotGap,subPlotH,subPlotW)
plot(output1.CGIters,'-','color',colorMatrix(1,:));
box on; grid on; hold on
plot(output2.CGIters,'color',colorMatrix(2,:));
plot(outputStar.CGIters,'color',colorMatrix(3,:));
xlabel('LDIPM Iteration','interpreter','Latex','FontSize',labelsize)
ylabel('CG Iterations','interpreter','Latex','FontSize',labelsize)
xlim([1 Inf])

subtightplot(3,1,3,subPlotGap,subPlotH,subPlotW)
h2 = semilogy([1 length(output2.condVec(:,1))],[CGTol CGTol],'linewidth',1,'color',colorMatrix(1,:));
box on; grid on; hold on
h3 = semilogy(output2.condVec(:,1),'linewidth',1,'color',colorMatrix(2,:));
% h4 = semilogy(output2.condVec(:,2),'linewidth',1,'color',colorMatrix(3,:));
xlabel('LDIPM Iteration','interpreter','Latex','FontSize',labelsize)
ylabel('Residual','interpreter','Latex','FontSize',labelsize)
legend('Nominal Requirement','Inexact Requirement','interpreter','Latex','FontSize',legendsize,'location','northeast')
% legend('Base Requirement','Requirement 1','Requirement 2','interpreter','Latex','FontSize',legendsize,'location','southwest')
ylim([1e-9 Inf])
xlim([1 Inf])
switch caseFlag
    case 3
        yticks([1e-8 1e-4 1e0])
    case 4
        yticks([1e-8 1e-4 1e0])
    case 6
        yticks([1e-8 1e-4 1e0])
end




if saveFigFlag == 1
    filename = strcat('./Plots/comparison',num2str(caseFlag),'_TRIPresPlot');
    saveas(gcf,filename,'epsc')
end
