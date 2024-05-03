clc
clear all
close all

% Load data
saveFlag = 0;
n = 1000;
m = 2000;
condNum = 100;
NSample = 1;


% LDIPM settings
mu_f = 1e-4;
mu_0 = 1e6;
maxIter = 150;
maxCGIter = 15000;
CGTol = 1e-6;
v0 = zeros(m,1);

% Set options
opts.mu_f = mu_f;
opts.mu_0 = mu_0;
opts.maxIter = maxIter;
opts.printFlag = 0;
opts.maxCGIter = maxCGIter;
opts.CGTol = CGTol;
opts.epsilon = 1e-3;

% Set parameters
gamma = 0.9;
opts.gamma = gamma;
opts.kappa = 0.1;
a_rand = 0.9;
b_rand = 1.1;

%% Main sampling iteration loop

% rng("default"); % set seed to startup seed

% Initialize stored data
CG_ITERS = zeros(NSample,2);
LDIPM_ITERS = zeros(NSample,2);
FEASFLAG = zeros(NSample,1);

% Run
for i = 1:NSample
    i
    % ////////////////////////////////////////
    % * Generate the QP matrices
    % ////////////////////////////////////////
    A = zeros(m,n);
    for ii = 1:m
        p_i = normrnd(zeros(1,n),1);
        a_i = p_i/norm(p_i);
        A(ii,:) = a_i;
    end

    % Get hh as a random NxN matrix
    r = n;
    R = zeros(r,n);
    for ii = 1:r
        p_i = normrnd(zeros(1,n),1);
        r_i = p_i/norm(p_i);
        R(ii,:) = r_i;
    end
    HPrime = R'*R;
    [u, s, v] = svd(HPrime);
    s = diag(s);           % s is vector
    s = s(1)*( 1-((condNum-1)/condNum)*(s(1)-s)/(s(1)-s(end))) ;
    s = diag(s);           % back to matrix
    H = u*s*v';
    H = 1/2*(H' + H);

    % Get the vectors
    x = normrnd(zeros(n,1),1);
    w1 = normrnd(zeros(m,1),1);
    w2 = normrnd(zeros(m,1),1);
%     s = ones(m,1) + 0.1*abs(w1);
    lambda = ones(m,1) + 0.1*abs(w2);
%     b = s - A*x;
    b = ones(m,1) + 0.1*abs(w1);
    while min(b) < 0
        b = ones(m,1) + 0.1*abs(w1);
    end
    c = A'*lambda - H*x;

    % ////////////////////////////////////////
    % * Get centered point for initialization
    % ////////////////////////////////////////
    dTol =  1e-4;
    [x_c,v_c,d] = logInteriorPoint_getCenteredPoint(H,c,A,b,v0,mu_0,1000,dTol);
    r_init =  a_rand + (b_rand-a_rand).*rand(m,1);
%     v_init = v_c.*r_init*0;
    v_init = v_c;

    % ////////////////////////////////////////
    % * Run both LDIPM algorithms
    % ////////////////////////////////////////
    [x1,output1] = cgLDIPM_longstep(H,c,A,b,v_init,opts);
    [x2,output2] = cgLDIPM_longstep_inexact(H,c,A,b,v_init,opts);

    % Process and store
    CG_ITERS(i,:) = [sum(sum(output1.CGIters)) sum(sum(output2.CGIters))];
    LDIPM_ITERS(i,:) = [output1.numIter output2.numIter];
    FEASFLAG(i) = output2.feasFlag;

%         % Save the first data set for a closer look
%     if i == 1  && saveFlag
%         saveStr = ['./Data/longstepSample_n',num2str(n),'_m',num2str(m),'_cond',num2str(condNum)];
%         sampleData.H = H;
%         sampleData.c = c;
%         sampleData.A = A;
%         sampleData.b = b;
%         sampleData.v_init = v_init;
%         sampleData.opts = opts;
%         sampleData.output1 = output1;
%         sampleData.output2 = output2;
%         save(saveStr,'sampleData');
%     end
end

% Save the data
if saveFlag
    saveData.CG_ITERS = CG_ITERS;
    saveData.LDIPM_ITERS = LDIPM_ITERS;
    saveData.n = n;
    saveData.m = m;
    saveData.FEASFLAG = FEASFLAG;
    saveData.condNum = condNum;
    saveData.NSample = NSample;
    saveStr = ['./Data/longstepData_n',num2str(n),'_m',num2str(m),'_cond',num2str(condNum)];
    save(saveStr,'saveData');
end

%% Run single processing step
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
legend([h1 h2],'Exact','Inexact','interpreter','Latex','FontSize',legendsize)
yticks([1e-4 1e-2 1e0 1e2 1e4])
xlim([1 Inf])
ylim([0.5*opts.mu_f Inf])
 
subtightplot(3,1,2,subPlotGap,subPlotH,subPlotW)
semilogy(output1.CGIters(:,1)+output1.CGIters(:,2),'color',colorMatrix(1,:));
box on; grid on; hold on
semilogy(output2.CGIters(:,1)+output2.CGIters(:,2),'color',colorMatrix(2,:),'linestyle','-.');
% legend('Exact','Inexact','interpreter','Latex','FontSize',legendsize,'location','southeast')
xlabel('$\mu$','interpreter','Latex','FontSize',labelsize)
ylabel('CG Iterations','interpreter','Latex','FontSize',labelsize)
xlim([1 Inf])
yticks([1e0 1e1 1e2 1e3 1e4])
ylim([5e1 2e4])

subtightplot(3,1,3,subPlotGap,subPlotH,subPlotW)
semilogy(ones(output1.numIter,1)*opts.CGTol,'-','color',colorMatrix(1,:));
box on; grid on; hold on
semilogy(output2.CGres(:,2),'color',colorMatrix(2,:),'linestyle',':');
semilogy(output2.CGres(:,1),'color',colorMatrix(2,:),'linestyle','-.');
legend('Exact','Inexact (1)','Inexact (2)','interpreter','Latex','FontSize',legendsize,'location','southwest','orientation','horizontal')
xlabel('LDIPM Iteration','interpreter','Latex','FontSize',labelsize)
ylabel('$\|r\|$','interpreter','Latex','FontSize',labelsize)
xlim([1 Inf])
ylim([0.25e-6 1e2])
yticks([1e-6 1e-4 1e-2 1e0 1e2])

