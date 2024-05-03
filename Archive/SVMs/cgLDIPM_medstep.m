function [x,output] = cgLDIPM_medStep(W,c,A,b,v0,opts)
% min 0.5*x'*W*x + c'*x   subject to:  A*x <= b

% Get size variables
m = size(A,1);
n = size(A,2);

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

if isfield(opts,'maxIter')
    maxIter = opts.maxIter;
else
    maxIter = 100;
end
% First, change variables to Ax + b >= 0... This is just for uniformity
% with quadprog's inputs.
A = - A;

% Pack
numIter = 0; % number of newton iterations performed

% --------------------- MAIN ITERATION LOOP ---------------------
% Unpack shortstep parameters
v = v0;
x = zeros(n,1);
mu = mu_0;
kappa  =  opts.kappa;
gamma = opts.gamma;

% Store CG output feedback
CGIters = zeros(maxIter,1);
CGres = zeros(maxIter,1);
muVec = zeros(maxIter,1);
count = 0;
while mu > mu_f
   mu  = (1/kappa)*mu;
   if mu < mu_f
       mu = mu_f;
   end
   dNorm = 10;
   while dNorm > 1
       [x,d,cgIters_i,res_i,] =  solveNewtonStep(mu,v,W,c,A,b,maxCGIter,CGTol,x);
       t = min([1, 1/(norm(d,'inf')^2)]);
       v = v + t*d;
       count = count + 1;
       CGIters(count) = cgIters_i;
       CGres(count)  = res_i;
       muVec(count) = mu;
       dNorm = norm(d,'inf');
   end
end

% Set output
output.v = v;
output.CGIters = CGIters(1:count,:);
output.CGres = CGres(1:count,:);
output.muVec = muVec(1:count);
output.numIter = count;
output.maxIter = maxIter;
end



function [x,d,numIter,res] = solveNewtonStep(mu,v,W,c,A,b,maxIter,tol,x0)
m = size(A,1);

% Define the RHS vector b
Q = diag(exp(2*v));
M = A'*Q*A + W;
f = 2*sqrt(mu)*A'*exp(v) - (c + A'*Q*b);
% mDiag = diag(M);

% --------------- CONJUGATE GRADIENT ---------------
% Initialize and redefine the problem such that x0 = 0.

% Run the first iteration of CG and iniialize iteration variables
x = x0;
r = f - M*x0;

% Calculate iteration constants
z = r;
% z = r./mDiag;
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
while numIter  < maxIter && res > tol 
    zPrev = z;
    z = r;
    % z = r./mDiag;
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
end
d = ones(m,1) - 1/sqrt(mu)*exp(v).*(A*x + b);
end




