clc
clear all
close all
addpath('/Users/jordan/Google Drive/SCHOOL/PhD/MATLAB/OptimizationFunctions')
pathVar = pathdef;
addpath(pathVar);

% Load data
saveFlag = 0;
caseFlag = 3;
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
    case 3 % Case 1 but with only the input constraints
        load('QPData');
        H = H_QP;
        c = f_QP;
        xmin = -ones(size(H_QP,1),1);
        xmax = ones(size(H_QP,1),1);
        [A,b] = minmaxMatrices(xmin,xmax);
    case 4
        N = 10;
        H = diag(linspace(0.01,1,N));
        %         H = eye(N);
        load('minmax')
        xmax = xmax*1000;
        xmin = xmin*1000;
        [A,b] = minmaxMatrices(xmin,xmax);        
end
invH = inv(H);


% Run LDIPM with normal settings
mu_f = 1e-14;
mu_0 = 1e8;
maxIter = 100;
maxCGIter = 10000;
CGPreCondFlag = 0;
printFlag = 0;
v0 = zeros(size(A,1),1);

xTolVec = logspace(-8,-2); 
xTolVec = fliplr(xTolVec);
numSample = length(xTolVec);
execTime_reg = zeros(numSample,1);
execTime_cg = zeros(numSample,1);
execTime_cgdiag = zeros(numSample,1);
execTime_pg = zeros(numSample,1);
numIters_reg = zeros(numSample,1);
numIters_cg = zeros(numSample,1);
numIters_cgdiag = zeros(numSample,1);
numIters_pg = zeros(numSample,1);
NumAverage = 200;


% Regular
brokeFlag = 0;
for iOuter = 1:numSample
    % Get optimal solution 
    xStar = logInteriorPoint(H,c,A,b,[],[],mu_f,mu_0,v0,maxIter,printFlag);
    xTol = xTolVec(iOuter);
    
    % Run and average
    execTimeSum = 0;
    for j = 1:NumAverage
        [~,xError,execTime_i,numIter] = logInteriorPoint_rt(H,c,A,b,mu_f,mu_0,v0,maxIter,xStar,xTol);
        execTimeSum = execTimeSum + execTime_i;
        if xError > xTol
            brokeFlag = 1;
            break
        end
    end
    if brokeFlag == 1
        fprintf('Error unachievable for Regular at xTol: %0.2e \n',xTol)
        fprintf('Produced Error: %0.2e \n',xError)
        numSample = iOuter-1;
        break
    else
        execTime_reg(iOuter) = execTimeSum/NumAverage;
        numIters_reg(iOuter) = numIter;
    end
end


% CG
brokeFlag = 0;
for iOuter = 1:numSample
    % Get optimal solution 
    xStar = logInteriorPoint(H,c,A,b,[],[],mu_f,mu_0,v0,maxIter,printFlag);
    xTol = xTolVec(iOuter);
    
    % Run and average
    execTimeSum = 0;
    for j = 1:NumAverage
        [~,xError,execTime_i,numIter] = logInteriorPoint_conjgrad_rt(H,c,A,b,mu_f,mu_0,v0,maxIter,maxCGIter,0,xStar,xTol);
        execTimeSum = execTimeSum + execTime_i;
        if xError > xTol
            brokeFlag = 1;
            break
        end
    end
    if brokeFlag == 1
        fprintf('Error unachievable for CG at xTol: %0.2e \n',xTol)
        fprintf('Produced Error: %0.2e \n',xError)
        numSample = iOuter-1;
        break
    else
        execTime_cg(iOuter) = execTimeSum/NumAverage;
        numIters_cg(iOuter) = numIter;
    end
end

% CG Diag
brokeFlag = 0;
for iOuter = 1:numSample
    % Get optimal solution 
    xStar = logInteriorPoint(H,c,A,b,[],[],mu_f,mu_0,v0,maxIter,printFlag);
    xTol = xTolVec(iOuter);
    
    % Run and average
    execTimeSum = 0;
    for j = 1:NumAverage
        [~,xError,execTime_i,numIter] = logInteriorPoint_conjgrad_rt(H,c,A,b,mu_f,mu_0,v0,maxIter,maxCGIter,1,xStar,xTol);
        execTimeSum = execTimeSum + execTime_i;
        if xError > xTol
            brokeFlag = 1;
            break
        end
    end
    if brokeFlag == 1
        fprintf('Error unachievable for CGdiag at xTol: %0.2e \n',xTol)
        fprintf('Produced Error: %0.2e \n',xError)
        numSample = iOuter-1;
        break
    else
        execTime_cgdiag(iOuter) = execTimeSum/NumAverage;
        numIters_cgdiag(iOuter) = numIter;
    end
end


% PG
if caseFlag == 3
    % Proj grad
    brokeFlag = 0;
    for iOuter = 1:numSample
        % Get optimal solution
        xStar = logInteriorPoint(H,c,A,b,[],[],mu_f,mu_0,v0,maxIter,printFlag);
        xTol = xTolVec(iOuter);
        
        % Run and average
        execTimeSum = 0;
        for j = 1:NumAverage
            [xx,numIter,xError,execTime_i] = projGradSolver_rt(H,c,zeros(size(H,1),1),xmin,xmax,50000,xStar,xTol);
            execTimeSum = execTimeSum + execTime_i;
            if xError > xTol
                brokeFlag = 1;
                break
            end
        end
        if brokeFlag == 1
            fprintf('Error unachievable for Proj Grad at xTol: %0.2e \n',xTol)
            fprintf('Produced Error: %0.2e \n',xError)
            numSample = iOuter-1;
            break
        else
            execTime_pg(iOuter) = execTimeSum/NumAverage;
            numIters_pg(iOuter) = numIter;
        end
    end
end

% Cut
xTolVec = xTolVec(1:numSample);
execTime_reg = execTime_reg(1:numSample);
execTime_cg = execTime_cg(1:numSample);
execTime_cgdiag = execTime_cgdiag(1:numSample);
numIters_reg = numIters_reg(1:numSample);
numIters_cg = numIters_cg(1:numSample);
numIters_cgdiag = numIters_cgdiag(1:numSample);
if caseFlag == 3
    execTime_pg = execTime_pg(1:numSample);
    numIters_pg = numIters_pg(1:numSample);
end


% dataname = ['.\Data\muPlot_Case',num2str(caseFlag)];
% save(dataname,'mu_f_vec','execTime_reg','execTime_cg')

%% Plotting
close all

figure
h1 = loglog(xTolVec,execTime_reg,'linewidth',2);
hold on; grid on; box on
h2 = loglog(xTolVec,execTime_cg,'linewidth',2);
h3 = loglog(xTolVec,execTime_cgdiag,'linewidth',2);
h4 = loglog(xTolVec,execTime_pg,'linewidth',2);
figSize = [0 0 0.3 0.3];
set(gcf,'units','normalized','position',figSize)
xlabel('$\| x - x^* \|$','interpreter','latex','fontsize',15)
ylabel('Exec Time (s)','interpreter','latex','fontsize',15)
legend([h1 h2 h3 h4],'Direct','CG','CG Diag','Grad','interpreter','latex','fontsize',12,'location','best')
% filename = [pwd '/Figures/xErrorPlot_Case',num2str(caseFlag)];
% saveas(gcf,filename,'epsc')

% figure
% ind = find(mu_f_vec <
% h1 = semilogx(1./mu_f_vec,execTime_reg,'linewidth',2);
% hold on; grid on; box on
% h2 = semilogx(1./mu_f_vec,execTime_cg,'linewidth',2);
% figSize = [0 0 0.2 0.2];
% set(gcf,'units','normalized','position',figSize)
% xlabel('$1/\mu_f$','interpreter','latex','fontsize',15)
% ylabel('Exec Time','interpreter','latex','fontsize',15)
% legend([h1 h2],'Reg','Diag','interpreter','latex','fontsize',12,'location','best')
% ylim([0 Inf])
