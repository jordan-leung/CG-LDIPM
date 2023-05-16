clc
clear all
close all

% Load data
saveFlag = 0;
caseFlag = 5;
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
mu_f = 1e-4;
mu_0 = 1;
maxIter = 150;
maxCGIter = 5;
CGTol = 0.1;
v0 = zeros(size(A,1),1);

% Set options
opts.mu_f = mu_f;
opts.mu_0 = mu_0;
opts.maxIter = maxIter;
opts.printFlag = 0;
opts.maxCGIter = maxCGIter;
opts.CGTol = CGTol;


% Get centered point
fprintf('------ Acquiring centred point  ------ \n')
opts_init = opts;
opts_init.mu_f = mu_0;
[x_init,output_init] = logInteriorPoint(H,c,A,b,v0,opts);
v0 = output_init.v;
[x_c,v_c,d] = logInteriorPoint_getCenteredPoint(H,c,A,b,v0,mu_0,1000,1e-8);
if norm(d,'inf') > 1e-6
    error('Didnt get initial centred point')
end
opts.x_c = x_c;
opts.v_c = v_c;

% Get initial condition somewhat close to v_c
% b = 0.2;
% a = -0.2;
% myRandVec = a + (b-a).*rand(length(v_c),1);
% save('myRandVec','myRandVec');
load('myRandVec');
v0 = v_c + 5*myRandVec.*v_c;

% % Run analysis
% fprintf('------ Running analysis  ------ \n')
% [x,output] = logInteriorPoint_conjgrad_analysis(H,c,A,b,v0,opts);

% Run analysis
fprintf('------ Running analysis  ------ \n')
[x,output] = logInteriorPoint_conjgrad_manual(H,c,A,b,v0,opts);

%% Plotting
close all
set(0,'defaultAxesFontSize',12)
set(0,'defaultLineLineWidth',2)
fontSize = 15;
legendSize = 12;
figSize = [0 0 0.3 0.3];

numberOfNans = sum(isnan(output.v_to_vr))

figure
set(gcf,'units','normalized','position',figSize)
subplot(2,2,1)
semilogy(output.v_to_vc)
hold on; box on; grid on
semilogy(output.v_to_vr,'-.')
semilogy(output.vr_to_vc,'--')
legend('$\| v - v_c \|$','$\| v - v_r \|$','$\| v_r - v_c \|$','interpreter','latex','fontsize',legendSize)
ylabel('$v$ Distances','interpreter','latex','fontsize',fontSize)

subplot(2,2,2)
semilogy(output.x_to_xc)
hold on; box on; grid on
semilogy(output.x_to_xr,'-.')
semilogy(output.xr_to_xc,'--')
legend('$\| x - x_c \|$','$\| x - x_r \|$','$\| x_r - x_c \|$','interpreter','latex','fontsize',legendSize)
ylabel('$x$ Distances','interpreter','latex','fontsize',fontSize)

subplot(2,2,3)
plot(output.dNorm)
hold on; box on; grid on
plot(output.dTrueNorm,'-.')
legend('$d$','$d^*$','interpreter','latex','interpreter','latex','fontsize',legendSize)
ylabel('$\| d\|$','interpreter','latex','fontsize',fontSize)

subplot(2,2,4)
semilogy(output.dError)
hold on; box on; grid on
semilogy(output.resNorm,'-.')
legend('$\| d - d^* \|$','$\| r \|$','interpreter','latex','fontsize',legendSize)
ylabel('CG Error','interpreter','latex','fontsize',fontSize)


saveFigFlag = 0;
if saveFigFlag
    filename = strcat('./Figures/','analysis_constRes');
    saveas(gcf,filename,'epsc'); 
end









