clc
clear all
close all
addpath('/Users/jordan/Google Drive/SCHOOL/PhD/MATLAB/OptimizationFunctions')
pathVar = pathdef;
addpath(pathVar);

% Load data
saveFlag = 0;
caseFlag = 1;
preCondFlag = 3;
switch caseFlag
    case 1
        load('QPData');
        H = H_QP;
        c = f_QP;
        A = A_QP;
        b = b_QP;
    case 2
        load('QPData2');
        H = H_QP;
        c = f_QP;
        A = A_QP;
        b = b_QP;       
end

% OPTIONS = optimoptions('quadprog');
% OPTIONS = optimoptions(OPTIONS, 'OptimalityTolerance', 1e-10, 'ConstraintTolerance', 1e-10);
% x_QP = quadprog(H,c,A,b,[],[],[],[],[],OPTIONS)

% Run LDIPM with normal settings
mu_f = 1e-10;
mu_0 = 1e8;
v0 = zeros(length(b),1);
% v0 = -10*ones(length(b),1);
maxIter = 1000;
maxCGIter = 1000; 
printFlag = 0;
[xStar,~,~,vStar,muStar,~,numIterStar,~,~,execTimeStar,CGIters,CGres] = logInteriorPoint_conjgrad(H,c,A,b,mu_f,mu_0,v0,maxIter,maxCGIter,preCondFlag,printFlag);


% Iterate
CGIterVec = 2:15;
p = length(CGIterVec);
numIterVec = zeros(p,1);
xError = zeros(p,1);
vError = zeros(p,1);
constraintSlack = zeros(p,1);
iter = 1;
for i = CGIterVec
    [x,~,~,v,~,~,numIter,~,~,execTimeStar,CGIters_i,CGres_i] = logInteriorPoint_conjgrad(H,c,A,b,mu_f,mu_0,v0,maxIter,i,preCondFlag,printFlag);
    
    numIterVec(iter) = numIter;
    xError(iter) = norm(x - xStar);
    vError(iter) = norm(v - vStar);
    constraintSlack(iter) = min(b - A*x);
    iter = iter + 1;
end
badConstraints = find(constraintSlack < 0);


%% Plotting
close all

figSize = [0 0 0.2 0.2];

% Errors
figure
set(gcf,'units','normalized','position',figSize)
h1 = semilogy(CGIterVec,xError,'.','Markersize',15);
hold on; grid on; box on;
semilogy(CGIterVec,xError,'Color',h1.Color)
h2 = semilogy(CGIterVec,vError,'.','Markersize',15);
semilogy(CGIterVec,vError,'Color',h2.Color)
xlabel('\# of CG Iterations','interpreter','latex','fontsize',15)
ylabel('Error','interpreter','latex','fontsize',15)
xlim([1 Inf])
legend([h1 h2],'$\|x - x^*\|$','$\|v - v^*\|$','interpreter','latex','fontsize',15)
if saveFlag
    switch caseFlag
        case 1
            switch preCondFlag
                case 0
                    saveas(gcf,'error_active','epsc')
                case 1
                    saveas(gcf,'error_active_preCond','epsc')
                case 3
                    saveas(gcf,'error_active_preCond_block','epsc')
            end
        case 2
            switch preCondFlag
                case 0
                    saveas(gcf,'error_inactive','epsc')
                case 1
                    saveas(gcf,'error_inactive_preCond','epsc')
                case 3
                    saveas(gcf,'error_inactive_preCond_block','epsc')
            end
    end
end

% Number of LDIPM Iterations
figure
set(gcf,'units','normalized','position',figSize)
h3 = plot(CGIterVec,numIterVec,'.','Markersize',15);
hold on; grid on; box on
plot(CGIterVec,numIterVec,'Color',h3.Color)
xlabel('\# of CG Iterations','interpreter','latex','fontsize',15)
ylabel('\# of LDIPM Iterations','interpreter','latex','fontsize',15)
xlim([1 Inf])
if saveFlag
    switch caseFlag
        case 1
            switch preCondFlag
                case 0
                    saveas(gcf,'iter_active','epsc')
                case 1
                    saveas(gcf,'iter_active_preCond','epsc')
                case 3
                    saveas(gcf,'iter_active_preCond_block','epsc')
            end
        case 2
            yticks([17 18])
            switch preCondFlag
                case 0
                    saveas(gcf,'iter_inactive','epsc')
                case 1
                    saveas(gcf,'iter_inactive_preCond','epsc')
                case 3
                    saveas(gcf,'iter_inactive_preCond_block','epsc')
            end
    end
end

% Constraint Slack
figure
set(gcf,'units','normalized','position',figSize)
h4 = plot(CGIterVec,constraintSlack,'.','Markersize',15);
hold on; grid on; box on
plot(CGIterVec,constraintSlack,'Color',h4.Color)
xlabel('\# of CG Iterations','interpreter','latex','fontsize',15)
ylabel('$\min_i b - Ax$','interpreter','latex','fontsize',15)
xlim([1 Inf])
if saveFlag
    switch caseFlag
        case 1
            switch preCondFlag
                case 0
                    saveas(gcf,'slack_active','epsc')
                case 1
                    saveas(gcf,'slack_active_preCond','epsc')
                case 3
                    saveas(gcf,'slack_active_preCond_block','epsc')
            end
        case 2
            switch preCondFlag
                case 0
                    saveas(gcf,'slack_inactive','epsc')
                case 1
                    saveas(gcf,'slack_inactive_preCond','epsc')
                case 3
                    saveas(gcf,'slack_inactive_preCond_block','epsc')
            end
    end
end

% Number of CG Iterations needed to achieve 1e-12 tolerance
figure
set(gcf,'units','normalized','position',figSize)
subplot(2,1,1)
h5 = plot(CGIters,'.','Markersize',15);
hold on; grid on; box on
plot(CGIters,'color',h5.Color);
ylabel('CG iterations','interpreter','latex','fontsize',15)
subplot(2,1,2)
h6 = semilogy(CGres,'.','Markersize',15);
hold on; grid on; box on
semilogy(CGres,'color',h6.Color);
xlabel('LDIPM Iteration','interpreter','latex','fontsize',15)
ylabel('CG Residual','interpreter','latex','fontsize',15)
if saveFlag
    switch caseFlag
        case 1
            switch preCondFlag
                case 0
                    saveas(gcf,'CGIters_active','epsc')
                case 1
                    saveas(gcf,'CGIters_active_preCond','epsc')
                case 3
                    saveas(gcf,'CGIters_active_preCond_block','epsc')
            end
        case 2
            switch preCondFlag
                case 0
                    saveas(gcf,'CGIters_inactive','epsc')
                case 1
                    saveas(gcf,'CGIters_inactive_preCond','epsc')
                case 3
                    saveas(gcf,'CGIters_inactive_preCond_block','epsc')
            end
    end
end
