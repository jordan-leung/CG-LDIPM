clc
clear all
close all
addpath('/Users/jordan/Google Drive/SCHOOL/PhD/MATLAB/OptimizationFunctions')
pathVar = pathdef;
addpath(pathVar);

% Load data
saveFlag = 0;
caseFlag = 1;
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
    case 3
        N = 10;
        H = diag(linspace(0.01,1,N));
        %         H = eye(N);
        load('minmax')
        xmax = xmax*1000;
        xmin = xmin*1000;
        [A,b] = minmaxMatrices(xmin,xmax);        
end
invH = inv(H);


OPTIONS = optimoptions('quadprog');
OPTIONS = optimoptions(OPTIONS, 'OptimalityTolerance', 1e-10, 'ConstraintTolerance', 1e-10);
x_QP = quadprog(H,c,A,b,[],[],[],[],[],OPTIONS)

% % Remove any zero rows
% i = 1;
% while i <= size(A,1)
%     if sum(abs(A(i,:))) < 1e-10
%         A(i,:) = [];
%         b(i) = [];
%     else
%         i = i + 1;
%     end
% end
% 
% % Calculate what invH-norm of each row looks like
% valVec = zeros(size(A,1),1);
% for i = 1:size(A,1)
%     a_i = A(i,:)';
%     val = a_i'*inv(H)*a_i;
%     valVec(i) = val;
% end
% targetVal = 1;
% 
% % Rescale the constraint matrix so that the w-norm is targetVal
% APrev = A;
% bPrev = b;
% GPrev = APrev*invH*APrev';
% lambdaVec = zeros(size(A,1),1);
% for i = 1:size(A,1)
%     lambda = sqrt(targetVal/valVec(i));
%     A(i,:) = lambda*A(i,:);
%     b(i,:) = lambda*b(i);
%     lambdaVec(i) = lambda;
% end
% 
% G = A*invH*A';
% [V,D] = eig(G);
% V = V';
% dd = diag(D);
% dd_alt = 0.001*ones(length(dd),1);
% for i = 1:length(dd_alt)
%     if dd(i) > 1e-6
%         dd_alt(i) = 1./sqrt(dd(i));
%     end
% end
% P = V'*diag(dd_alt)*V;
% Pnew = diag(diag(P));
% Gnew = Pnew*A*invH*A'*Pnew;
% APrev = A;
% bPrev = b;
% A = Pnew*A;
% b = Pnew*b;


% Run LDIPM with normal settings
mu_f = 1e-10;
mu_0 = 1e8;
% v0 = (log(1e-3) - log(sqrt(mu_0)))*ones(size(A,1),1);
maxIter = 1000;
maxCGIter = 10000;
CGPreCondFlag = 0;
printFlag = 1;
v0 = zeros(size(A,1),1);

% % % No rescaling
% [x_reg,~,~,vStar_reg,muStar_reg,~,numIter_reg,~,~,~,CGIters_reg,CGres_reg] = logInteriorPoint_conjgrad(H,c,APrev,bPrev,mu_f,mu_0,v0,maxIter,maxCGIter,0,printFlag);

% No conditioning
fprintf('------ Running with default ------ \n')
[x,~,~,vStar,muStar,~,numIter,~,~,~,CGIters,CGres] = logInteriorPoint_conjgrad(H,c,A,b,mu_f,mu_0,v0,maxIter,maxCGIter,0,1);

fprintf('------ Running with simple scheme ------ \n')

% % Approximate newton
% [x1,~,~,v1,mu1,~,numIter1,~,~,~,CGIters1,CGres1,dError1] = logInteriorPoint_inexactNewton(H,c,A,b,mu_f,mu_0,v0,maxIter,maxCGIter,0,printFlag);

% % Diagonal
% [x2,~,~,v2,~,~,numIter2,~,~,~,CGIters2,CGres2] = logInteriorPoint_conjgrad(H,c,A,b,mu_f,mu_0,v0,maxIter,maxCGIter,1,printFlag);
% 
% % % Diagonal + rank 1
% % [x3,~,~,v3,~,~,numIter3,~,~,~,CGIters3,CGres3] = logInteriorPoint_conjgrad(H,c,A,b,mu_f,mu_0,v0,maxIter,maxCGIter,2,printFlag);
% 
% Block diagonal
[x4,~,~,v4,~,~,numIter4,~,~,~,CGIters4,CGres4] = logInteriorPoint_conjgrad(H,c,A,b,mu_f,mu_0,v0,maxIter,maxCGIter,3,printFlag);

[CGIters CGIters4]

% figure
% hold on; grid on; box on
% h1 = plot(CGIters,'.','Markersize',15);
% h2 = plot(CGIters2,'.','Markersize',15);
% % h3 = plot(CGIters3,'.','Markersize',15);
% h4 = plot(CGIters4,'.','Markersize',15);
% plot(CGIters,'Color',h1.Color)
% plot(CGIters2,'Color',h2.Color)
% % plot(CGIters3,'Color',h3.Color)
% plot(CGIters4,'Color',h4.Color)
% figSize = [0 0 0.2 0.2];
% set(gcf,'units','normalized','position',figSize)
% xlabel('LDIPM Iteration','interpreter','latex','fontsize',15)
% ylabel('\# of CG Iterations','interpreter','latex','fontsize',15)
% legend([h1 h2 h4],'None','Diag','Block','interpreter','latex','fontsize',12,'location','best')
% % legend([h1 h2],'None','Diag','interpreter','latex','fontsize',12,'location','best')
% ylim([0 Inf])


% [xReg,lambdaReg,sReg,vReg,~,~,numIterReg] = logInteriorPoint(H,c,A,b,[],[],mu_f,mu_0,v0,maxIter,printFlag);