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

% Fixed settings
CGTol = 1e-10;
preCondFlag = 1;
v0 = zeros(size(A,1),1);
if caseFlag < 5
    vThresh = -4;
    vNumThresh = length(v0)/4;
else
    vThresh = 0;
    vNumThresh = -1; % so always on
end

% Varying settings
mu = 1e-2;
numSamples = 10000;
[~,~,~,v0] = logInteriorPoint(H,c,A,b,[],[],mu,1e8,zeros(size(A,1),1),100000,0); % run to get v0

% Run
[fHist,dHist,fPrimeHist,fPrimeInitHist,CGIters,CGIters_thresh,CGIters_decrease] = divergenceAnalysis(H,c,A,b,mu,v0,CGTol,preCondFlag,vThresh,vNumThresh,numSamples);

%% Plotting
close all
muString = num2str(mu,'%.E');

% Only plot stuff less than 10
indVec_plot = find(fHist < 10);

figure
set(gcf,'units','normalized','position',[0.1300 0.1192 0.1 0.2])
fPrimeBound = -fHist - dHist.^2;
plot(fHist(indVec_plot),fPrimeHist(indVec_plot),'.')
hold on; box on; grid on
plot(fHist(indVec_plot),fPrimeInitHist(indVec_plot),'.')
plot(fHist(indVec_plot),fPrimeBound(indVec_plot),'.')
plot(xlim,[0 0],'k')
xlabel('$f(0)$','interpreter','latex','fontsize',15)
ylabel('$f''(0)$','interpreter','latex','fontsize',15)
legend('f''(0)','f''(0) after 1 CG iteration','$-f(0)-\|d\|^2$','interpreter','latex','fontsize',12,'location','best')
% saveas(gcf,['divTraj_',muString],'png')


% Find points where fPrimeInitHit > 0
indVec_pos = find(fPrimeInitHist > 0);
indVec = intersect(indVec_pos,indVec_plot);


% Get the values of f at this
fPos = fHist(indVec);

% Generate the same histogram 
figure
set(gcf,'units','normalized','position',[0.1300 0.1192 0.1 0.2])
subplot(2,1,1)
h1 = histogram(fHist(indVec_plot));
grid on; box on;
xlabel('$f(0)$','interpreter','latex','fontsize',15)
ylabel('Total Samples','interpreter','latex','fontsize',15)
subplot(2,1,2)
h2 = histogram(fPos,h1.BinEdges);
grid on; box on;
xlabel('$f(0)$','interpreter','latex','fontsize',15)
ylabel('$f''(0) > 0$','interpreter','latex','fontsize',15)
% saveas(gcf,['divHistogram_mu_',muString],'png')

figure
plot(fHist(indVec),fPrimeInitHist(indVec),'.')
