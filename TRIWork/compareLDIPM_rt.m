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
maxIter = 1;
maxCGIter = 1;
CGPreCondFlag = 0;
printFlag = 0;
v0 = zeros(size(A,1),1);
xTol = 1e-6;


%% Compile

useMexFlag = 1;
compileFlag = 1;
if compileFlag
    fprintf('----------- Compiling code ---------- \n')
    if caseFlag == 3
        codegen projGradSolver_rt -args {H,c,ones(size(H,1),1),xmin,xmax,50000,ones(size(H,1),1),1}
    end
    codegen logInteriorPoint_rt -args {H,c,A,b,mu_f,mu_0,v0,maxIter,ones(size(H,1),1),1}
    codegen logInteriorPoint_conjgrad_rt -args {H,c,A,b,mu_f,mu_0,v0,maxIter,maxCGIter,0,ones(size(H,1),1),1}
end


% Get optimal solution
xStar = logInteriorPoint(H,c,A,b,[],[],mu_f,mu_0,v0,200,printFlag);

%% Run

% PG
if caseFlag == 3
    if useMexFlag
        [x_pg,numIter_pg,xError_pg,execTime_pg] = projGradSolver_rt_mex(H,c,zeros(size(H,1),1),xmin,xmax,50000,xStar,xTol);
    else
        [x_pg,numIter_pg,xError_pg,execTime_pg] = projGradSolver_rt(H,c,zeros(size(H,1),1),xmin,xmax,50000,xStar,xTol);
    end
end

% Regular
if useMexFlag
    [x_reg,xError_reg,execTime_reg,numIter_reg] = logInteriorPoint_rt_mex(H,c,A,b,mu_f,mu_0,v0,maxIter,xStar,xTol);
else
    [x_reg,xError_reg,execTime_reg,numIter_reg] = logInteriorPoint_rt(H,c,A,b,mu_f,mu_0,v0,maxIter,xStar,xTol);
end

% CG
if useMexFlag
    [x_cg,xError_cg,execTime_cg,numIter_cg] = logInteriorPoint_conjgrad_rt_mex(H,c,A,b,mu_f,mu_0,v0,maxIter,maxCGIter,0,xStar,xTol);
else
    [x_cg,xError_cg,execTime_cg,numIter_cg] = logInteriorPoint_conjgrad_rt(H,c,A,b,mu_f,mu_0,v0,maxIter,maxCGIter,0,xStar,xTol);
end

% CG
if useMexFlag
    [~,xError_cg_full,execTime_cg_full,numIter_cg_full] = logInteriorPoint_conjgrad_rt_mex(H,c,A,b,mu_f,mu_0,v0,maxIter,10000,0,xStar,xTol);
else
    [~,xError_cg_full,execTime_cg_full,numIter_cg_full] = logInteriorPoint_conjgrad_rt(H,c,A,b,mu_f,mu_0,v0,maxIter,10000,0,xStar,xTol);
end
% % CG Diag
% [x_diag,xError_diag,execTime_diag,numIter_diag] = logInteriorPoint_conjgrad_rt(H,c,A,b,mu_f,mu_0,v0,maxIter,maxCGIter,1,xStar,xTol);


time = [execTime_reg execTime_cg execTime_cg_full]

%% Plot error
close all


figure
semilogy(xError_cg);
hold on; grid on; box on
semilogy(xError_cg_full);
xlabel('LDIPM Iteration','interpreter','latex')
ylabel('$\|x - x^*\|$','interpreter','latex')
legend('MaxCG = 6','No max')

figure
semilogy(xError_reg)


