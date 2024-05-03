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
n = 140;
m = n/2;

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
N = 10;
R = 0.1*eye(m);
nonzeros = round(.7*n);
dTol = 1e-3;

% Storage containers
CG_ITERS = zeros(NSample,2);
LDIPM_ITERS = zeros(NSample,2);
FEASFLAG = zeros(NSample,1);

% Saveflag
saveFlag = 0;

for i = 1:NSample
    i
    % (A,B) matrices
    A = eye(n) + normrnd(zeros(n,n),0.01*ones(n,n));
    A = 0.99*A/max(abs(eig(A)));
    B = normrnd(zeros(n,m),1*ones(n,m));

    % (Q,R) matrices
    q = zeros(n,1);
    q(1:nonzeros) = 10*rand(nonzeros,1);
    q = q(randperm(n));
    Q = diag(q);
    [K,P] = dlqr(A,B,Q,R);

    % Constraints
    aa = 1*ones(n,1);
    bb = 2*ones(n,1);
    xBar =  aa + (bb-aa).*rand(n,1);
    uBar = 0.1*rand(m,1);

    % Initial Condition
    aa = -0.5*xBar;
    bb = 0.5*xBar;
    x0 =  aa + (bb-aa).*rand(n,1);

    % MPC and QPmatrices 
    [W,F,Acon,L,b] = generateMatrices_forCGLDIPM(N,A,B,P,Q,R,xBar,uBar);
    if i == 1
        size(Acon)
    end
    c = 2*F*x0;
    bcon = b + L*x0;
    v0 = zeros(size(Acon,1),1);

    % Get centered point 
    [~,v_c,d,iters] = logInteriorPoint_getCenteredPoint(W,c,Acon,bcon,v0,opts.mu_0,100,dTol);
    if iters == 100
        warning('Did not arrive at centered point')
    end    
    % Run exact CG
    [x1,output1] = cgLDIPM_longstep(W,c,Acon,bcon,v_c,opts);
    % Run inexact CG
    [x2,output2] = cgLDIPM_longstep_inexact(W,c,Acon,bcon,v_c,opts);

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

iLast = NSample;
% iLast = i-1;
% Save the data
if saveFlag
    saveData.CG_ITERS = CG_ITERS(1:iLast,:);
    saveData.LDIPM_ITERS = LDIPM_ITERS(1:iLast,:);
    saveData.AVG_ITERS = CG_ITERS(1:iLast,:)./LDIPM_ITERS(1:iLast,:);
    saveData.n = n;
    saveData.m = m;
    saveData.n_qp = size(Acon,2);
    saveData.m_qp = size(Acon,1);
    saveData.FEASFLAG = FEASFLAG(1:iLast);
    saveData.NSample = iLast;
    saveStr = ['./Data/LongStep/longstepData_n',num2str(n)];
    save(saveStr,'saveData');
else
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
    semilogy(length(muVecTotal1),muVecTotal1(end),'.' ,'color',colorMatrix(1,:),'markersize',18);
    h2 = semilogy(muVecTotal2,'color',colorMatrix(2,:),'linestyle','-.');
    semilogy(length(muVecTotal2),muVecTotal2(end),'.' ,'color',colorMatrix(2,:),'markersize',18);
    xlabel('Cumulative CG Iterations','interpreter','Latex','FontSize',labelsize)
    ylabel('$\mu$','interpreter','Latex','FontSize',labelsize)
    legend([h1 h2],'Nominal','Inexact','interpreter','Latex','FontSize',legendsize)
    yticks([1e-4 1e0 1e4])
    xlim([1 Inf])

    subtightplot(3,1,2,subPlotGap,subPlotH,subPlotW)
    semilogy(output1.CGIters(:,1)+output1.CGIters(:,2),'color',colorMatrix(1,:));
    box on; grid on; hold on
    semilogy(output2.CGIters(:,1)+output2.CGIters(:,2),'color',colorMatrix(2,:),'linestyle','-.');
    legend('Nominal','Inexact','interpreter','Latex','FontSize',legendsize,'location','southeast')
    xlabel('$\mu$','interpreter','Latex','FontSize',labelsize)
    ylabel('CG Iterations','interpreter','Latex','FontSize',labelsize)
    xlim([1 Inf])
    yticks([1e0 1e2 1e4])

    subtightplot(3,1,3,subPlotGap,subPlotH,subPlotW)
    semilogy(ones(output1.numIter,1)*opts.CGTol,'-','color',colorMatrix(1,:));
    box on; grid on; hold on
    semilogy(output2.CGres(:,2),'color',colorMatrix(2,:),'linestyle',':');
    semilogy(output2.CGres(:,1),'color',colorMatrix(2,:),'linestyle','-.');
    legend('Nominal','Inexact (1)','Inexact (2)','interpreter','Latex','FontSize',legendsize,'location','southwest','orientation','horizontal')
    xlabel('LDIPM Iteration','interpreter','Latex','FontSize',labelsize)
    ylabel('Required Residual','interpreter','Latex','FontSize',labelsize)
    xlim([1 Inf])
    ylim([0.25e-6 Inf])
    yticks([1e-6 1e-3 1e0 1e3])
end




