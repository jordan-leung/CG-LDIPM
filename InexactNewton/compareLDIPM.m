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
        H = H_QP;CG
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
mu_f = 1e-7;
mu_0 = 1e8;
maxIter = 50000;
maxCGIter = 100000;
CGTol = 1e-8;
v0 = zeros(size(A,1),1);

% Set options
opts.mu_f = mu_f;
opts.mu_0 = mu_0;
opts.maxIter = maxIter;
opts.printFlag = 1;
opts.maxCGIter = maxCGIter;
opts.CGTol = CGTol;
opts.preCondFlag = 0;


fprintf('------ Running with CG LDIPM (longstep)  ------ \n')
[x_cg,output_cg] = logInteriorPoint_conjgrad(H,c,A,b,v0,opts);

% fprintf('------ Running with CGLDIPM-GIN  ------ \n')
% opts.FNormLimit = 5;
% [x1,output1] = logInteriorPoint_conjgrad_INB(H,c,A,b,v0,opts);
% 
% fprintf('------ Running with CGLDIPM-GIN  ------ \n')
% opts.FNormLimit = 1;
% [x2,output2] = logInteriorPoint_conjgrad_INB(H,c,A,b,v0,opts);
% 
% fprintf('------ Running with CGLDIPM-GIN  ------ \n')
% opts.FNormLimit = 0.5;
% [x3,output3] = logInteriorPoint_conjgrad_INB(H,c,A,b,v0,opts);

fprintf('------ Running with CGLDIPM-GIN  ------ \n')
opts.FNormLimit = 5;
[x1,output1] = logInteriorPoint_conjgrad_INB_edit(H,c,A,b,v0,opts);

fprintf('------ Running with CGLDIPM-GIN  ------ \n')
opts.FNormLimit = 1;
[x2,output2] = logInteriorPoint_conjgrad_INB_edit(H,c,A,b,v0,opts);

fprintf('------ Running with CGLDIPM-GIN  ------ \n')
opts.FNormLimit = 0.5;
[x3,output3] = logInteriorPoint_conjgrad_INB_edit(H,c,A,b,v0,opts);


% Solution check
norms = zeros(3,1);
norms(1) = norm(x_cg - x1);
norms(2)= norm(x_cg - x2);
norms(3) = norm(x_cg - x3);
for i = 1:length(norms)
   if norms(i) > 1e-4 || isnan(norms(i))
       warning('Big norm!!!!')
       fprintf('Solution %0.0d has a norm of %0.2e \n',i,norms(i))
   end
end


%% PLotting
close all
set(0,'defaultAxesFontSize',12)
figSize = [0 0 0.2 0.2];
saveFigFlag = 0;

% Mu vs. CG iterations 
figure
totalIter_cg = cumsum((sum(output_cg.CGIters'))');
totalIter1 = cumsum(output1.CGIters);
totalIter2 = cumsum(output2.CGIters);
totalIter3 = cumsum(output3.CGIters);
set(gcf,'units','normalized','position',figSize)
h1 = semilogy(totalIter_cg,output_cg.muVec);
hold on
h2 = semilogy(totalIter1,output1.muVec);
h3 = semilogy(totalIter2,output2.muVec);
h4 = semilogy(totalIter3,output3.muVec);
semilogy(totalIter_cg,output_cg.muVec,'.','Markersize',15,'Color',h1.Color)
semilogy(totalIter1,output1.muVec,'.','Markersize',15,'Color',h2.Color)
semilogy(totalIter2,output2.muVec,'.','Markersize',15,'Color',h3.Color)
semilogy(totalIter3,output3.muVec,'.','Markersize',15,'Color',h4.Color)
semilogy(totalIter_cg(end),output_cg.muVec(end),'s','Markersize',15,'Linewidth',2,'Color',h1.Color)
semilogy(totalIter1(end),output1.muVec(end),'s','Markersize',15,'Linewidth',2,'Color',h2.Color)
semilogy(totalIter2(end),output1.muVec(end),'s','Markersize',15,'Linewidth',2,'Color',h3.Color)
semilogy(totalIter3(end),output1.muVec(end),'s','Markersize',15,'Linewidth',2,'Color',h4.Color)
grid on; box on;
hold off
legend([h1 h2 h3 h4],'Longstep w/ CG','Inexact Newton, $\epsilon = 10$',...
    'Inexact Newton, $\epsilon = 1$','Inexact Newton, $\epsilon = 0.5$',...
    'interpreter','latex','fontsize',12,'location','northeast')
xlabel('CG Iterations','interpreter','latex','fontsize',15)
ylabel('$\mu$','interpreter','latex','fontsize',15)
if saveFigFlag
    filename = strcat('./Figures/','preliminaryPlot');
    saveas(gcf,filename,'epsc'); 
end

% % CG iterations vs. LDIPM iteration
% set(0, 'DefaultLineLineWidth', 2);
% figure
% CGIters_longstep =  output_cg.CGIters(:,1) + output_cg.CGIters(:,2);
% set(gcf,'units','normalized','position',figSize)
% plot(CGIters_longstep,'color',h1.Color);
% grid on; box on; hold on
% plot(output1.CGIters,'color',h2.Color);
% plot(output2.CGIters,'color',h3.Color);
% plot(output3.CGIters,'color',h4.Color);
% plot(length(CGIters_longstep),CGIters_longstep(end),'s','Markersize',15,'Linewidth',2,'Color',h1.Color)
% plot(length(output1.CGIters),output1.CGIters(end),'s','Markersize',15,'Linewidth',2,'Color',h2.Color)
% plot(length(output2.CGIters),output2.CGIters(end),'s','Markersize',15,'Linewidth',2,'Color',h3.Color)
% plot(length(output3.CGIters),output3.CGIters(end),'s','Markersize',15,'Linewidth',2,'Color',h4.Color)
% legend('Longstep w/ CG','Inexact Newton, $\epsilon = 10$',...
%     'Inexact Newton, $\epsilon = 1$','Inexact Newton, $\epsilon = 0.5$',...
%     'interpreter','latex','fontsize',12,'location','northeast')
% xlabel('LDIPM Iterations','interpreter','latex','fontsize',15)
% ylabel('CG Iteratios','interpreter','latex','fontsize',15)
% if saveFigFlag
%     filename = strcat('./Figures/','preliminaryPlot_iters');
%     saveas(gcf,filename,'epsc'); 
% end
% 
% % Feasibility
% figure
% set(gcf,'units','normalized','position',figSize)
% plot(output1.feasVec,'color',h2.Color);
% hold on; box on; grid on;
% plot(output2.feasVec,'color',h3.Color);
% plot(output3.feasVec,'color',h4.Color);
% plot(length(output1.feasVec),output1.feasVec(end),'s','Markersize',15,'Linewidth',2,'Color',h2.Color)
% plot(length(output2.feasVec),output2.feasVec(end),'s','Markersize',15,'Linewidth',2,'Color',h3.Color)
% plot(length(output3.feasVec),output3.feasVec(end),'s','Markersize',15,'Linewidth',2,'Color',h4.Color)
% xlabel('LDIPM Iterations','interpreter','latex','fontsize',15)
% ylabel('Primal Feasibility','interpreter','latex','fontsize',15)
% if saveFigFlag
%     filename = strcat('./Figures/','preliminaryPlot_feas');
%     saveas(gcf,filename,'epsc'); 
% end
% 

