clc
clear all
close all

% Load data
saveFlag = 0;
caseFlag = 1;
switch caseFlag
    case 1 % Active constraints MPC example
        load('QPData'); 
        H = H_QP;
        c = f_QP;
        A = A_QP;
        b = b_QP;
    case 2 % Inactive constraints MPC example
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
    case 4 % Scalable QP example
        N = 10;
        H = diag(linspace(0.01,1,N));
        %         H = eye(N);
        load('minmax')
        [A,b] = minmaxMatrices(xmin,xmax);
    case 5 % Much more constraints than variables
        load('Case5_Data');
end
invH = inv(H);

% LDIPM settings
mu_f = 1e-10;
mu_0 = 1e8;
maxIter = 50000;
maxCGIter = 100000;
CGTol = 1e-10;
v0 = zeros(size(A,1),1);

% Shortstep_optimal - no WS
fprintf('------ Running with CGLDIPM-GIN  ------ \n')
[x,lambda,s,v,mu,startFlag,numIter,muStar,exitFlag,execTime,CGIters,CGres] = logInteriorPoint_conjgrad_INB(H,c,A,b,mu_f,mu_0,v0,maxIter,maxCGIter,CGTol,1);


