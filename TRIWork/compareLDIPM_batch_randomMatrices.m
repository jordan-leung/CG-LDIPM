clc
clear all
close all
addpath('/Users/jordan/Google Drive/SCHOOL/PhD/MATLAB/OptimizationFunctions')
pathVar = pathdef;
addpath(pathVar);

% Load data - mainly for compiling code... Will regenerate H, c, b in the
% main for loop

% Static variables
N = 120;
condTarget = 1e3;
ratio = 1;
m = N/ratio;
A = zeros(m,N);
for i = 1:m
    A(i,1+(i-1)*(N/m):i*(N/m)) = ones(1,N/m);
end
b_low = -1;
b_high = 1;
c_low = -10;
c_high = 10;

% Changing Variables
h1= 2*rand(N,N)-1 + 2*rand(N,N)-1;
h2 = h1*h1';              % symmetric with random entries beween -2 and 2
[u, s, v] = svd(h2);
s = diag(s);           % s is vector
s = s(1)*( 1-((condTarget-1)/condTarget)*(s(1)-s)/(s(1)-s(end))) ;
s = diag(s);           % back to matrix
H = u*s*v';
H = 1/2*(H' + H); 
c =  c_low + (c_high-c_low).*rand(N,1);
b = b_low + (b_high - b_low)*rand(m,1);



% Run LDIPM with normal settings
mu_f = 1e-10;
mu_0 = 1e8;
maxIter = 100;
maxCGIter = 1000;
CGTol = 1e-2;
CGPreCondFlag = 0;
printFlag = 0;
v0 = zeros(size(A,1),1);

% Define sample space
numSample = 100;
execTime_reg = zeros(numSample,1);
execTime_cg = zeros(numSample,1);
CGIters = zeros(numSample,1);
IPMIters  = zeros(numSample,1);
NumAverage = 25;

%% Produce the codegen files

useMexFlag = 1;
compileFlag = 1;
if compileFlag && useMexFlag
    fprintf('----------- Compiling code ---------- \n')
    codegen logInteriorPoint_rt -args {H,c,A,b,mu_f,mu_0,v0,maxIter}
    codegen logInteriorPoint_conjgrad_rt -args {H,c,A,b,mu_f,mu_0,v0,maxIter,CGTol,maxCGIter,0}
end

%% Run experiment


% Regular
fprintf('----------- Running regular ---------- \n')
brokeFlag = 0;
iOuter = 1;
while iOuter <= numSample
    if mod(iOuter,10) == 0
    fprintf('Iteration number: %0.0i \n',iOuter);
    end
    
    % Generate problem
    hh= 2*rand(N,N)-1 + 2*rand(N,N)-1;
    hh = hh*hh';              % symmetric with random entries beween -2 and 2
    [u, s, v] = svd(hh);
    s = diag(s);           % s is vector
    s = s(1)*( 1-((condTarget-1)/condTarget)*(s(1)-s)/(s(1)-s(end))) ;
    s = diag(s);           % back to matrix
    H = u*s*v';
    H = 1/2*(H' + H);
    c =  c_low + (c_high-c_low).*rand(N,1);
    b = b_low + (b_high - b_low)*rand(m,1);
    
    % Get optimal solution
    OPTIONS = optimoptions('quadprog');
    OPTIONS = optimoptions(OPTIONS, 'OptimalityTolerance', 1e-12, 'ConstraintTolerance', 1e-10,'Display','off');
    xStar = quadprog(H,c,A,b,[],[],[],[],[],OPTIONS);
    
    % Run and average
    execTimeSum_reg = 0;
    execTimeSum_cg = 0;
    workFlag = 1;
    for j = 1:NumAverage
        if useMexFlag
            [x1,mu1,execTime_1,numIter1] = logInteriorPoint_rt_mex(H,c,A,b,mu_f,mu_0,v0,maxIter);
            [x2,mu2,execTime_2,numIter2,CGIters_i] = logInteriorPoint_conjgrad_rt_mex(H,c,A,b,mu_f,mu_0,v0,maxIter,maxCGIter,CGTol,1);
        else
            [x1,mu1,execTime_1,numIter1] = logInteriorPoint_rt(H,c,A,b,mu_f,mu_0,v0,maxIter);
            [x2,mu2,execTime_2,numIter2,CGIters_i] = logInteriorPoint_conjgrad_rt(H,c,A,b,mu_f,mu_0,v0,maxIter,maxCGIter,CGTol,1);
        end
        if norm(x1 - xStar) < 1e-6 && norm(x2 - xStar) < 1e-6
            execTimeSum_reg = execTimeSum_reg + execTime_1;
            execTimeSum_cg = execTimeSum_cg + execTime_2;
        else
            fprintf('Broken at iteration %0.0i \n',iOuter);
            fprtinf('x1Error = %0.2e, x2Error = %0.2e \n',norm(x1 - xStar),norm(x2 - xStar) );
            workFlag = 0;
            break
        end
    end
    if workFlag
        execTime_reg(iOuter) = execTimeSum_reg/NumAverage;
        execTime_cg(iOuter) = execTimeSum_cg/NumAverage;
        CGIters(iOuter) = CGIters_i;
        IPMIters(iOuter) = numIter2;
        iOuter = iOuter + 1;
    end % otherwise we repeat
end



% dataname = ['.\Data\muPlot_Case',num2str(caseFlag)];
% save(dataname,'mu_f_vec','execTime_reg','execTime_cg')

%% Plotting
close all
saveFlag = 1;

figure
semilogy(execTime_reg,'r.')
hold on; box on; grid on
semilogy(execTime_cg,'b.')
legend('LDIPM','CG-LDIPM','location','Best')
ylabel('Execution Time')
xlabel('Trial')
figSize = [0 0 0.2 0.2];
set(gcf,'units','normalized','position',figSize)
set(gcf,'Resize','off')
if saveFlag
    saveas(gcf,['./Figures/randomMatrices_exec_',num2str(ratio)],'epsc')
end

figure
plot(IPMIters,'b.')
ylabel('IPM Iterations')
xlabel('Trial')
grid on; box on
figSize = [0 0 0.2 0.2];
set(gcf,'units','normalized','position',figSize)
set(gcf,'Resize','off')
if saveFlag
    saveas(gcf,['./Figures/randomMatrices_IPM_',num2str(ratio)],'epsc')
end

figure
plot(CGIters,'b.')
ylabel('Total CG Iterations')
xlabel('Trial')
grid on; box on
figSize = [0 0 0.2 0.2];
set(gcf,'units','normalized','position',figSize)
set(gcf,'Resize','off')
if saveFlag
    saveas(gcf,['./Figures/randomMatrices_CG_',num2str(ratio)],'epsc')
end

figure
plot(CGIters./IPMIters,'b.')
ylabel('Average CG Iterations')
xlabel('Trial')
grid on; box on
figSize = [0 0 0.2 0.2];
set(gcf,'units','normalized','position',figSize)
set(gcf,'Resize','off')
if saveFlag
    saveas(gcf,['./Figures/randomMatrices_CGNorm_',num2str(ratio)],'epsc')
end

