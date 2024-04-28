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
mu_f = 1e-2;
mu_0 = 1e4;
maxIter = 1e4;
maxCGIter = 1000;
CGTol = 1e-6;
v0 = zeros(size(A,1),1);

% Set options
opts.mu_f = mu_f;
opts.mu_0 = mu_0;
opts.maxIter = maxIter;
opts.printFlag = 0;
opts.maxCGIter = maxCGIter;
opts.CGTol = CGTol;
opts.beta = 0.45;

fprintf('------ Running analysis  ------ \n')
[x_ss,output_ss] = logInteriorPoint_conjgrad_shortstep(H,c,A,b,v0,opts);

fprintf('------ Running analysis  ------ \n')
[x,output] = logInteriorPoint_conjgrad_cheat(H,c,A,b,v0,opts);

norm(x-x_ss)

% Post-process to get plot of mu (y-axis) vs. total cg iterations
CGIters = output.CGIters;
muVec = output.muVec;
totalIters = sum(CGIters);
muVec_total = zeros(totalIters,1);
count = 1;
for i = 1:length(CGIters)
    for j = 1:CGIters(i)
        muVec_total(count) = muVec(i);
        count = count + 1;
    end
end

% Repeat for regular  ss
CGIters_ss = output_ss.CGIters;
muVec_ss = output_ss.muVec;
totalIters_ss = sum(CGIters_ss);
muVec_total_ss = zeros(totalIters_ss,1);
count = 1;
for i = 1:length(CGIters_ss)
    for j = 1:CGIters_ss(i)
        muVec_total_ss(count) = muVec_ss(i);
        count = count + 1;
    end
end


%% Plotting
close all
set(0,'defaultAxesFontSize',12)
set(0,'defaultLineLineWidth',2)
fontSize = 15;
legendSize = 12;
figSize = [0 0 0.3 0.3];


figure
set(gcf,'units','normalized','position',figSize)
plot(output.CGIters)
hold on; box on; grid on
plot(output_ss.CGIters)
legend('Inexact CG','Exact CG','interpreter','latex','fontsize',legendSize)
xlabel('LDIPM Iterations','interpreter','latex','fontsize',legendSize)
ylabel('CG Iterations','interpreter','latex','fontsize',fontSize)

figure
set(gcf,'units','normalized','position',figSize)
semilogy(muVec_total)
hold on; box on; grid on
semilogy(muVec_total_ss)
legend('Inexact CG','Exact CG','interpreter','latex','fontsize',legendSize)
xlabel('CG Iterations','interpreter','latex','fontsize',fontSize)
ylabel('$\mu$','interpreter','latex','fontsize',legendSize)

figure
set(gcf,'units','normalized','position',figSize)
plot(output.h_vHat_v)
hold on; box on; grid on
plot(output_ss.h_vHat_v)
legend('Inexact CG','Exact CG','interpreter','latex','fontsize',legendSize)
xlabel('LDIPM Iterations','interpreter','latex','fontsize',legendSize)
ylabel('$h_\mu(v)$','interpreter','latex','fontsize',fontSize)


figure
set(gcf,'units','normalized','position',figSize)
semilogy(output.resNorm)
hold on; box on; grid on
semilogy(output_ss.resNorm)
legend('Inexact CG','Exact CG','interpreter','latex','fontsize',legendSize)
xlabel('LDIPM Iterations','interpreter','latex','fontsize',legendSize)
ylabel('$\|r\|$','interpreter','latex','fontsize',fontSize)

