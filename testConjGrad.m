clc
clear all
close all

N = 500;

% % Generate random parameters and save
% P = RandOrthMat(N);
% aa = 1;
% bb = 0.1;
% WVec = aa + (bb-aa).*rand(N-1,1);
% aa = 5;
% bb = -5;
% b =  aa + (bb-aa).*rand(N,1);
% % save('data','P','WVec','b');


% % Run
% % load('data')
% W1 = diag([1000; WVec]);
% W2 = (P\W1)*P;
% W = W2*W2';
% [U,S,V] = svd(W);
% u = U(:,1);
% v = V(:,1);
% GMatrix = S(1,1)*u*v';
% G = [S(1,1)*u v];
% A = eye(N) + W;

% % Generate random parameters and save
% load('data')
% aa = 5;
% bb = -5;
% w =  aa + (bb-aa).*rand(N,1);
% GMatrix = w*w';
% G = [w w];
% A = eye(N) + GMatrix;



x = A\b;
[x1,numIter1,resStar1] = conjGrad(A,b,x0,maxIter,tol);
     
       
% [x2,numIter2,resStar2] = conjGrad_preCond(A,b,x0,200,tol,2,G);
% 
% [x3,numIter3,resStar3] = conjGrad_preCond(A,b,x0,200,tol,1);


numIter1
numIter2
numIter3