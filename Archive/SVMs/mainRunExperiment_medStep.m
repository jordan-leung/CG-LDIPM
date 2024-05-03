clear
close all
clc

addpath('/Users/jordan/Documents/GitHub/optimizationFunctions');
addpath('/Users/jordan/Documents/GitHub/extraFunctions')
addpath('/Users/jordan/Documents/GitHub/CG-LDIPM/MPC_Examples')

set(groot,'defaultFigureCreateFcn',@(fig,~)addToolbarExplorationButtons(fig))

%% * Set simulation parameters

% Set dimensions
% targetProblemSize = 100;
% n = 2*targetProblemSize/10
n = 100;
m = 2*n;

% Set number of trials to average over
NSample = 1;

% LDIPM settings
opts.mu_f =  1e-4;
opts.mu_0 = 1e6;
opts.maxIter = 100;
opts.printFlag = 0;
opts.maxCGIter = 20000;
opts.CGTol = 1e-6;
opts.gamma = 0.9;
opts.N = 1;
opts.kappa = (1/0.2);

% Parameters that are constant between trials
dTol = 1e-3;
nonzeros = round(0.15*n);
W = eye(n);
c = 10*ones(n,1);

% Storage containers
CG_ITERS = zeros(NSample,2);
LDIPM_ITERS = zeros(NSample,2);
FEASFLAG = zeros(NSample,1);

% Saveflag
saveFlag = 0;

for i = 1:NSample
    i
    % Define the vector b
    b = [ones(m/2,1); -1*ones(m/2,1)];
    A = zeros(m,n);
    % Fill in 15% non-zero elements
    for count = 1:m/2
      a_i = normrnd(1/n*ones(n,1),1/n*ones(n,1));
      a_i(nonzeros+1:end) = 0*a_i(nonzeros+1:end);
      a_i = a_i(randperm(n));
      A(count,:) = a_i';
    end
    for count = m/2+1:m
      a_i = normrnd(-1/n*ones(n,1),1/n*ones(n,1));
      a_i(nonzeros+1:end) = 0*a_i(nonzeros+1:end);
      a_i = a_i(randperm(n));
      A(count,:) = a_i';
    end
    v0 = zeros(size(A,1),1);

    % Generate ACon and bCon
    bcon = [b; -b.*(A*)]

    % Get centered point
    [~,v_c,d,iters] = logInteriorPoint_getCenteredPoint(W,c,A,b,v0,opts.mu_0,100,dTol);
    if iters == 100
        warning('Did not arrive at centered point')
    end
    % Run exact CG
    [x1,output1] = cgLDIPM_medstep(W,c,A,b,v_c,opts);
    % Run inexact CG
    [x2,output2] = cgLDIPM_medstep_inexact(W,c,A,b,v_c,opts);

    % Process and store
    CG_ITERS(i,:) = [sum(sum(output1.CGIters)) sum(sum(output2.CGIters))];
    LDIPM_ITERS(i,:) = [output1.numIter output2.numIter];
    FEASFLAG(i) = output2.feasFlag;

    %         % Save the first data set for a closer look
    %         if i == 1  && saveFlag
    %             saveStr = ['./Data/fixedstepSample_n',num2str(n),'_m',num2str(m),'_cond',num2str(condNum)];
    %             sampleData.H = H;
    %             sampleData.c = c;
    %             sampleData.A = A;
    %             sampleData.b = b;
    %             sampleData.v_init = v_init;
    %             sampleData.opts = opts;
    %             sampleData.output1 = output1;
    %             sampleData.output2 = output2;
    %             save(saveStr,'sampleData');
    %         end

end

x_qp = quadprog(W,c,A,b);

iLast = NSample;
% iLast = i-1;
% Save the data
if saveFlag
    saveData.CG_ITERS = CG_ITERS(1:iLast,:);
    saveData.LDIPM_ITERS = LDIPM_ITERS(1:iLast,:);
    saveData.AVG_ITERS = CG_ITERS(1:iLast,:)./LDIPM_ITERS(1:iLast,:);
    saveData.n = n;
    saveData.m = m;
    saveData.n_qp = size(A,2);
    saveData.m_qp = size(A,1);
    saveData.FEASFLAG = FEASFLAG(1:iLast);
    saveData.NSample = iLast;
    saveStr = ['./Data/FixedStep/fixedstepData_n',num2str(n)];
    save(saveStr,'saveData');
else
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
    CGIters1 = CGIters;

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

    %% Plotting
    close all
    set(0,'defaultLineLineWidth',2)
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
    semilogy(length(muVecTotal1),muVecTotal1(end),'.','color',colorMatrix(1,:),'markersize',18);
    h2 = semilogy(muVecTotal2,'color',colorMatrix(2,:),'linestyle','-.');
    semilogy(length(muVecTotal2),muVecTotal2(end),'.' ,'color',colorMatrix(2,:),'markersize',18);
    xlabel('Cumulative CG Iterations','interpreter','Latex','FontSize',labelsize)
    ylabel('$\mu$','interpreter','Latex','FontSize',labelsize)
    legend([h1 h2],'Nominal','Inexact','interpreter','Latex','FontSize',legendsize)
    yticks([1e-4 1e0 1e4])
    xlim([1 Inf])

    subtightplot(3,1,2,subPlotGap,subPlotH,subPlotW)
    semilogy(output1.CGIters(:,1),'color',colorMatrix(1,:));
    box on; grid on; hold on
    semilogy(output2.CGIters(:,1),'color',colorMatrix(2,:),'linestyle','-.');
    legend('Nominal','Inexact','interpreter','Latex','FontSize',legendsize,'location','southeast')
    xlabel('$\mu$','interpreter','Latex','FontSize',labelsize)
    ylabel('CG Iterations','interpreter','Latex','FontSize',labelsize)
    xlim([1 Inf])
    yticks([1e0 1e2 1e4])


    subtightplot(3,1,3,subPlotGap,subPlotH,subPlotW)
    semilogy(ones(output1.numIter,1)*opts.CGTol,'-','color',colorMatrix(1,:));
    box on; grid on; hold on
    semilogy(output2.resLim(:,1),'color',colorMatrix(2,:),'linestyle','-.');
    legend('Nominal','Inexact','interpreter','Latex','FontSize',legendsize)
    xlabel('LDIPM Iteration','interpreter','Latex','FontSize',labelsize)
    ylabel('Required Residual','interpreter','Latex','FontSize',labelsize)
    xlim([1 Inf])
    ylim([0.25e-6 Inf])
    yticks([1e-6 1e-3 1e0 1e3])
end

