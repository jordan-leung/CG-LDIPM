clc
clear all
close all

% Load data
saveFigFlag = 0;
caseFlag = 4;
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
mu_f = 1e-6;
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
opts.gamma = gamma;

% Get centered point
fprintf('------ Acquiring centred point  ------ \n')
opts_init = opts;
opts_init.mu_f = mu_0;
dTol =  1e-4;
[x_c,v_c,d] = logInteriorPoint_getCenteredPoint(H,c,A,b,v0,mu_0,1000,dTol);
if norm(d,'inf') > dTol
    error('Didnt get initial centred point')
end

 
fprintf('------ Running longstep LDIPM  ------ \n')
[x1,output1] = cgLDIPM_longstep(H,c,A,b,v_c,opts);
fprintf('------ Running inexact longstep CG-LDIPM  ------ \n')
[x2,output2] = cgLDIPM_longstep_inexact_noRestrict(H,c,A,b,v_c,opts);
% fprintf('------ Running inexact longstep CG-LDIPM  ------ \n')
% [x2,output2] = cgLDIPM_longstep_inexact(H,c,A,b,v_c,opts);


norm(x2-x1,'inf')
norm(output2.v - output1.v,'inf')

%% Plotting

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
totalIters = sum(CGIters); % gives 1 x 2 vector of summed columns
muVecTotal2 = zeros(sum(totalIters),1);
count = 1;
for i = 1:size(CGIters,1)
    for j = 1:(CGIters(i,1)+CGIters(i,2))
        muVecTotal2(count) = output2.muVec(i);
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
xlabel('Cumulative CG Iterations','interpreter','Latex','FontSize',labelsize)
ylabel('$\mu$','interpreter','Latex','FontSize',labelsize)
legend([h1 h2],'Nominal','Inexact','interpreter','Latex','FontSize',legendsize)
yticks([1e-4 1e0 1e4])
xlim([1 Inf])

subtightplot(3,1,2,subPlotGap,subPlotH,subPlotW)
plot(output1.CGIters,'-','color',colorMatrix(1,:));
box on; grid on; hold on
plot(output2.CGIters(:,2),'color',colorMatrix(2,:),'linestyle',':');
plot(output2.CGIters(:,1),'color',colorMatrix(2,:),'linestyle','-.');
plot(output2.CGIters(:,1)+output2.CGIters(:,2),'color',colorMatrix(2,:),'linestyle','-');
legend('Nominal','Inexact (1)','Inexact (2)','Inexact (Total)','interpreter','Latex','FontSize',legendsize)
xlabel('LDIPM Iteration','interpreter','Latex','FontSize',labelsize)
ylabel('CG Iterations','interpreter','Latex','FontSize',labelsize)
xlim([1 Inf])

subtightplot(3,1,3,subPlotGap,subPlotH,subPlotW)
semilogy(ones(output1.numIter,1)*opts.CGTol,'-','color',colorMatrix(1,:));
box on; grid on; hold on
semilogy(output2.CGres(:,2),'color',colorMatrix(2,:),'linestyle',':');
semilogy(output2.CGres(:,1),'color',colorMatrix(2,:),'linestyle','-.');
legend('Nominal','Inexact (1)','Inexact (2)','interpreter','Latex','FontSize',legendsize)
xlabel('LDIPM Iteration','interpreter','Latex','FontSize',labelsize)
ylabel('Required Residual','interpreter','Latex','FontSize',labelsize)
xlim([1 Inf])





if saveFigFlag == 1
    filename = strcat('./Plots/comparison',num2str(caseFlag),'_longstep');
    saveas(gcf,filename,'epsc')
end
