function [x,output] = cgLDIPM_inexact(W,c,Aineq,bineq,v0,opts)
% min 0.5*x'*W*x + c'*x   subject to:  A*x <= b

% Get size variables
m = size(Aineq,1);
n = size(Aineq,2);

% Set options
if isfield(opts,'mu_f')
    mu_f = opts.mu_f;
else
    mu_f = 1e-8;
end
if isfield(opts,'mu_0')
    mu_0 = opts.mu_0;
else
    mu_0 = 1e8;
end
if isfield(opts,'printFlag')
    printFlag = opts.printFlag;
else
    printFlag = 0;
end
if isfield(opts,'maxCGIter')
    maxCGIter = opts.maxCGIter;
else
    maxCGIter = 1000;
end
if isfield(opts,'CGTol')
    CGTol = opts.CGTol;
else
    CGTol = 1e-8;
end

% First, change variables to Ax + b >= 0... This is just for uniformity
% with quadprog's inputs.
Aineq = - Aineq;

% Pack
invW = inv(W);
const.W = W;
const.invW = invW;
const.c = c;
const.A = Aineq;
const.b = bineq;
numIter = 0; % number of newton iterations performed

% --------------------- MAIN ITERATION LOOP ---------------------
% Unpack shortstep parameters
v = v0;
x = zeros(n,1);
mu = mu_0;
kappa  =  opts.kappa;
N = opts.N;
const.gamma = opts.gamma;
const.minEig = min(eig(W));

% Store CG output feedback
maxIter = ceil(N*log(mu_0/mu_f)*sqrt(m)/(2*qInv(opts.zeta^2))); % just for initialization
CGIters = zeros(maxIter,1);
CGres = zeros(maxIter,1);
muVec = zeros(maxIter,1);
minCond = zeros(maxIter,1);
minInd = zeros(maxIter,1);
condVec = zeros(maxIter,2);
sigma = zeros(maxIter,1);
count = 0;
dNorm = 2;
while mu > mu_f && dNorm > 1
   if mu > mu_f
       mu  = (1/kappa)*mu;
   end
   for i = 1:N
       [x,d,cgIters_i,res_i,minCond_i,minInd_i,condVec_i,sigma_i] =  solveNewtonStep(mu,v,const,maxIter,CGTol,x);
       v = v + d;
       count = count + 1;
       CGIters(count) = cgIters_i;
       CGres(count)  = res_i;
       muVec(count) = mu;
       minCond(count) = minCond_i;
       minInd(count) = minInd_i;
       sigma(count) = sigma_i;
       condVec(count,:) = condVec_i'/sigma_i;
   end
   dNorm = norm(d,'inf');
end

% Set output
output.v = v;
output.CGIters = CGIters(1:count,:);
output.CGres = CGres(1:count,:);
output.muVec = muVec(1:count);
output.minCond = minCond(1:count);
output.minInd = minInd(1:count);
output.sigma = sigma(1:count);
output.condVec = condVec(1:count,:);
output.reqRes =  minCond(1:count)./sigma(1:count);
output.numIter = count;
output.maxIter = maxIter;
end



function [x,d,numIter,res,minCond,minInd,condVec,sigma] = solveNewtonStep(mu,v,const,maxIter,tol,x0)
W = const.W;
% invW = const.invW;
c = const.c;
A = const.A;
b = const.b;
m = size(A,1);
n = size(A,2);
gamma = const.gamma;
minEig = const.minEig;

% Define the RHS vector b
Q = diag(exp(2*v));
M = A'*Q*A + W;
f = 2*sqrt(mu)*A'*exp(v) - (c + A'*Q*b);

% --------------- CONJUGATE GRADIENT ---------------
% Initialize and redefine the problem such that x0 = 0.

% Run the first iteration of CG and iniialize iteration variables
x = x0;
r = f - M*x0;

% Calculate iteration constants
z = r;
p = z;
w = M*p;
alpha = r'*z/(p'*w);

% Update x and r
x = x + alpha*p; 
rPrev = r;
r = r - alpha*w; % r1
res = norm(r,2);
numIter = 1;

% Iterate...
truncCond = 0;


while numIter  < maxIter && res > tol && truncCond ==  0
    zPrev = z;
    z = r;
    beta = r'*z/(rPrev'*zPrev);
    p = z + beta*p;
    w = M*p;
    alpha = r'*z/(p'*w);
    
    % Update x
    x = x + alpha*p;
    rPrev = r;
    r = r - alpha*w;
    res = norm(r,2);

    % i++
    numIter = numIter + 1;

    % Check truncation criteria
    s_i = A*x + b;
    val_i = 0;
    for j = 1:m
        val_i = val_i + A(j,:)'/s_i(j);
    end
    sigma = (2/(minEig*mu))*norm(W*x + c - mu*val_i,2);
    d = ones(m,1) - 1/sqrt(mu)*exp(v).*s_i;
    dNorm = norm(d,2);
    cond1 = (1-gamma)*dNorm^2;
    if sigma*res < cond1
        gamma_i = 1 - sigma*res/dNorm^2;
    else
        gamma_i = gamma;
    end
    hlb_i = hlb(dNorm,gamma_i);
    cond2 = 0.5*(1 - gamma_i*sqrt(4*(1-gamma_i) + dNorm^2*gamma_i^2))*hlb_i^2;
    condVec = [cond1; cond2];
    [minCond,minInd] = min(condVec);
    if minCond > 0 && sigma*res < minCond
        truncCond = 1;
    end
end
end

% function val = barrierGrad(x,A,b)
%    m = size(A,1);
%    n = size(A,2);
%    val = zeros(n,1);
%    for i = 1:m
%        a_i = A(i,:)';
%        b_i = b(i);
%        val = val - a_i/(a_i'*x+b_i);
%    end
% end


function val = hlb(dNorm,gamma)
val  = gamma^2*dNorm^2/(2 - gamma + sqrt(4*(1-gamma) +  dNorm^2*gamma^2));
end
