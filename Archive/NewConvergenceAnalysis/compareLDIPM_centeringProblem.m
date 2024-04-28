clc
clear all
close all

% Load data
saveFigFlag = 0;
caseFlag = 3;
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
mu_0 = 1;
maxIter = 150;
maxCGIter = 50000;
CGTol = 1e-8;

% Set options
opts.mu_f = mu_f;
opts.mu_0 = mu_0;
opts.maxIter = maxIter;
opts.printFlag = 0;
opts.maxCGIter = maxCGIter;
opts.CGTol = CGTol;


% Get centered point
fprintf('------ Acquiring centred point  ------ \n')
v0 = zeros(size(A,1),1);
opts_init = opts;
opts_init.mu_f = mu_0;
opts_init.mu_0 = 1e8;
[x_c,v_c,d,v_init] = logInteriorPoint_getCenteredPoint(H,c,A,b,v0,mu_0,1000,1e-6);
if norm(d,'inf') > 1e-6
    error('Didnt get initial centred point')
end

% Get centered point using full C
% ar = 0;
% br = 2;
% randVec = ar + (br-ar)*rand(length(v_c),1); 
[x0,v0] = logInteriorPoint_getCenteredPoint(H,c,A,b,v0,mu_0*1e4,1000,1e-6);
[x1,v1,output1] = getCenteredPoint_base2(H,c,A,b,x0,v0,mu_0,1000,CGTol,v_c);
gamma = 0.9;
[x2,v2,output2] = getCenteredPoint_inexact(H,c,A,b,x0,v0,mu_0,1000,CGTol,v_c,gamma);
% [x2,v2,output2] = getCenteredPoint_inexact(H,c,A,b,x_c,v_c,mu_0,1000,CGTol,v_c,gamma);


%% Plotting
close all

% Post-process to get plot of mu (y-axis) vs. total cg iterations
CGIters = output1.numIterVec;
hVec = output1.hVec;
totalIters = sum(CGIters);
hVecTotal1 = zeros(totalIters,1);
count = 1;
for i = 1:length(CGIters)
    for j = 1:CGIters(i)
        hVecTotal1(count) = hVec(i);
        count = count + 1;
    end
end

% Post-process to get plot of mu (y-axis) vs. total cg iterations
CGIters = output2.numIterVec;
hVec = output2.hVec;
totalIters = sum(CGIters);
hVecTotal2 = zeros(totalIters,1);
count = 1;
for i = 1:length(CGIters)
    for j = 1:CGIters(i)
        hVecTotal2(count) = hVec(i);
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
semilogy(hVecTotal1,'color',colorMatrix(1,:));
box on; grid on; hold on
semilogy(hVecTotal2,'color',colorMatrix(2,:));
xlabel('Cumulative CG Iterations','interpreter','Latex','FontSize',labelsize)
ylabel('$h_\mu(v)$','interpreter','Latex','FontSize',labelsize)
legend('Base','Inexact','interpreter','Latex','FontSize',legendsize)
yticks([1e-8 1e-4 1e0 1e4])
ylim([1e-8 Inf])

subtightplot(4,1,2,subPlotGap,subPlotH,subPlotW)
plot(output1.numIterVec,'-','color',colorMatrix(1,:));
box on; grid on; hold on
plot(output2.numIterVec,'color',colorMatrix(2,:));
xlabel('LDIPM Iteration','interpreter','Latex','FontSize',labelsize)
ylabel('CG Iterations','interpreter','Latex','FontSize',labelsize)

subtightplot(4,1,3,subPlotGap,subPlotH,subPlotW)
h2 = semilogy([1 length(output2.delta)],[CGTol CGTol],'linewidth',1,'color',colorMatrix(1,:));
box on; grid on; hold on
h3 = semilogy(output2.delta./output2.sigma,'linewidth',1,'color',colorMatrix(2,:));
xlabel('LDIPM Iteration','interpreter','Latex','FontSize',labelsize)
ylabel('Residual','interpreter','Latex','FontSize',labelsize)
legend('Base Requirement','Requirement','interpreter','Latex','FontSize',legendsize,'location','southwest')
ylim([1e-9 Inf])
switch caseFlag
    case 3
        yticks([1e-8  1e-4  1e0])
    case 4
        yticks([1e-8  1e-4  1e0])
    case 6
        yticks([1e-8  1e-4  1e0])
end

subtightplot(4,2,7,subPlotGap,subPlotH,subPlotW)
semilogy(output2.sigma,'color',colorMatrix(2,:))
box on; grid on; hold on
xlabel('LDIPM Iteration','interpreter','Latex','FontSize',labelsize)
ylabel('$\sigma_{\mu}(x)$','interpreter','Latex','FontSize',labelsize)
switch caseFlag
    case 3
        yticks([1e-6 1e-3  1e0  1e3])
    case 4
        yticks([1e-6 1e-3  1e0  1e3])
    case 6
        yticks([1e-6 1e-3  1e0  1e3])
end

subtightplot(4,2,8,subPlotGap,subPlotH,subPlotW)
semilogy(output2.delta,'color',colorMatrix(2,:))
box on; grid on; hold on
xlabel('LDIPM Iteration','interpreter','Latex','FontSize',labelsize)
ylabel('$\delta(\|d\|)$','interpreter','Latex','FontSize',labelsize)
switch caseFlag
    case 3
        yticks([1e-10 1e-5  1e0  1e5])
    case 4
        yticks([1e-10 1e-5  1e0  1e5])
    case 6
        yticks([1e-10 1e-5  1e0  1e5])
end



if saveFigFlag == 1
    filename = strcat('./Plots/centering',num2str(caseFlag));
    saveas(gcf,filename,'epsc')
end
