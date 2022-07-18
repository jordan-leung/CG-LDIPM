clc
clear all
close all


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


% Static variables
N = 500;
condTarget = 1e3;
m = N;
b_low = -1;
b_high = 1;

% Changing Variables
hh= 2*rand(N,N)-1 + 2*rand(N,N)-1;
hh = hh*hh';              % symmetric with random entries beween -2 and 2
[u, s, v] = svd(hh);
s = diag(s);           % s is vector
s = s(1)*( 1-((condTarget-1)/condTarget)*(s(1)-s)/(s(1)-s(end))) ;
s = diag(s);           % back to matrix
H = u*s*v';
H = 1/2*(H' + H);
b = b_low + (b_high - b_low)*rand(m,1);
x0 = b*0;
maxIter = 5*N;
tol = 1e-10;
[x1,numIter1,resStar1] = conjGrad(H,b,x0,maxIter,tol);
norm(x1 - H\b) 

% A = magic(5);
% b = rand(5,1);
% x0 = b*0;
% maxIter = 50;
% tol = 1e-10;
% [x1,numIter1,resStar1] = conjGrad(A,b,x0,maxIter,tol);
% norm(x1 - A\b) 

       
% [x2,numIter2,resStar2] = conjGrad_preCond(A,b,x0,200,tol,2,G);
% 
% [x3,numIter3,resStar3] = conjGrad_preCond(A,b,x0,200,tol,1);

% 
% numIter1
% numIter2
% numIter3