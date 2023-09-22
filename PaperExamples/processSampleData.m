clc
clear all
close all
set(groot,'defaultFigureCreateFcn',@(fig,~)addToolbarExplorationButtons(fig))

n = 1000;
m = 2000;
condNum = 100;
saveFigFlag = 1;

%% Process fixed step

loadStr = ['./Data/fixedstepSample_n',num2str(n),'_m',num2str(m),'_cond',num2str(condNum)];
load(loadStr);
output1 = sampleData.output1;
output2 = sampleData.output2;
opts = sampleData.opts;

% Post-process to get plot of mu (y-axis) vs. total cg iterations
CGIters = output1.CGIters; % for base
totalIters = sum((CGIters));
muVecTotal1 = zeros(totalIters,1);
count = 1;
for i = 1:size(CGIters,1)
    for j = 1:CGIters(i,1)
        muVecTotal1(count) = output1.muVec(i);
        count = count + 1;
    end
end

CGIters = output2.CGIters; % for inexact
totalIters = sum((CGIters)); % gives 1 x 2 vector of summed columns
muVecTotal2 = zeros(totalIters,1);
count = 1;
for i = 1:size(CGIters,1)
    for j = 1:CGIters(i,1)
        muVecTotal2(count) = output2.muVec(i);
        count = count + 1;
    end
end


% Plotting
set(0,'defaultLineLineWidth',1)
set(0,'defaultAxesFontSize',15)
labelsize = 18;
legendsize = 15;
colorMatrix = colororder;
colorMatrix(7,:) = [0 102 0]/255;
subPlotGap = [0.1 0.1];
subPlotH = [0.1 0.05];
subPlotW = [0.15 0.05];
figSize = [0 0 0.25 0.35];

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
legend([h1 h2],'Exact','Inexact','interpreter','Latex','FontSize',legendsize)
yticks([1e-4 1e-2 1e0 1e2 1e4 1e6])
xlim([1 Inf])
ylim([0.5e-4 Inf])

subtightplot(3,1,2,subPlotGap,subPlotH,subPlotW)
semilogy(output1.CGIters(:,1),'color',colorMatrix(1,:));
box on; grid on; hold on
semilogy(output2.CGIters(:,1),'color',colorMatrix(2,:),'linestyle','-');
% legend('Exact','Inexact','interpreter','Latex','FontSize',legendsize)
xlabel('LDIPM Iteration','interpreter','Latex','FontSize',labelsize)
ylabel('CG Iterations','interpreter','Latex','FontSize',labelsize)
xlim([1 Inf])
ylim([1e1 1e4])
yticks([1e1 1e2 1e3 1e4])

subtightplot(3,1,3,subPlotGap,subPlotH,subPlotW)
semilogy(ones(output1.numIter,1)*opts.CGTol,'-','color',colorMatrix(1,:));
box on; grid on; hold on
semilogy(output2.CGres(:,1),'color',colorMatrix(2,:),'linestyle','-');
% legend('Exact','Inexact','interpreter','Latex','FontSize',legendsize)
xlabel('LDIPM Iteration','interpreter','Latex','FontSize',labelsize)
ylabel('$\|r\|$ Restriction','interpreter','Latex','FontSize',labelsize)
xlim([1 Inf])
ylim([0.5e-6 Inf])
yticks([1e-6 1e-4 1e-2 1e0 1e2])


if saveFigFlag == 1
    filename = strcat('./Plots/comparison_fixedstep');
    saveas(gcf,filename,'epsc')
end

%% Process longstep


loadStr = ['./Data/longstepSample_n',num2str(n),'_m',num2str(m),'_cond',num2str(condNum)];
load(loadStr);
output1 = sampleData.output1;
output2 = sampleData.output2;
opts = sampleData.opts;


% Post-process to get plot of mu (y-axis) vs. total cg iterations
CGIters = output1.CGIters; % for base
totalIters = sum(sum(CGIters));
muVecTotal1 = zeros(totalIters,1);
count = 1;
for i = 1:size(CGIters,1)
    for j = 1:(CGIters(i,1)+CGIters(i,2))
        muVecTotal1(count) = output1.muVec(i);
        count = count + 1;
    end
end
CGIters1 = CGIters;

CGIters = output2.CGIters; % for inexact
totalIters = sum(sum(CGIters)); % gives 1 x 2 vector of summed columns
muVecTotal2 = zeros(totalIters,1);
count = 1;
for i = 1:size(CGIters,1)
    for j = 1:(CGIters(i,1)+CGIters(i,2))
        muVecTotal2(count) = output2.muVec(i);
        count = count + 1;
    end
end


% Plotting

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
legend([h1 h2],'Exact','Inexact','interpreter','Latex','FontSize',legendsize)
yticks([1e-4 1e-2 1e0 1e2 1e4 1e6])
xlim([1 Inf])
ylim([0.5e-4 Inf])

subtightplot(3,1,2,subPlotGap,subPlotH,subPlotW)
semilogy(output1.CGIters(:,1)+output1.CGIters(:,2),'color',colorMatrix(1,:));
box on; grid on; hold on
semilogy(output2.CGIters(:,1)+output2.CGIters(:,2),'color',colorMatrix(2,:),'linestyle','-');
xlabel('LDIPM Iteration','interpreter','Latex','FontSize',labelsize)
ylabel('CG Iterations','interpreter','Latex','FontSize',labelsize)
xlim([1 Inf])
yticks([1e1 1e2 1e3 1e4])
ylim([1e2 1e4])

subtightplot(3,1,3,subPlotGap,subPlotH,subPlotW)
semilogy(ones(output1.numIter,1)*opts.CGTol,'-','color',colorMatrix(1,:));
box on; grid on; hold on
semilogy(output2.CGres(:,1),'color',colorMatrix(2,:),'linestyle','-.');
semilogy(output2.CGres(:,2),'color',colorMatrix(2,:),'linestyle','--');
legend('Exact','Inexact (1)','Inexact (2)','interpreter','Latex','FontSize',legendsize,'numcolumns',3)
xlabel('LDIPM Iteration','interpreter','Latex','FontSize',labelsize)
ylabel('$\|r\|$ Restriction','interpreter','Latex','FontSize',labelsize)
xlim([1 Inf])
ylim([0.5e-6 4e2])
yticks([1e-6 1e-4 1e-2 1e0 1e2])

if saveFigFlag == 1
    filename = strcat('./Plots/comparison_longstep');
    saveas(gcf,filename,'epsc')
end

