function [x,output] = cgLDIPM_medStep_inexact(W,c,Aineq,bineq,v0,opts)
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
if isfield(opts,'x0')
    x0 = opts.x0;
else
    x0 = zeros(n,1);
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
x = x0;
xFeas = x0;
mu = mu_0;
kappa  =  opts.kappa;
gamma = opts.gamma;
const.gamma = gamma;
const.minEig = min(eig(W));

% Store CG output feedback
maxIter = opts.maxIter; % just for initialization
CGIters = zeros(maxIter,1);
CGres = zeros(maxIter,1);
muVec = zeros(maxIter,1);
resLim = zeros(maxIter,1);
count = 0;
while mu > mu_f
    mu  = (1/kappa)*mu;
    if mu < mu_f
        mu = mu_f;
    end
    dNorm = 10;
    while dNorm > 1
        [x,d,cgIters_i,r_i,resLim_i,xFeas] =  solveNewtonStep(mu,v,const,maxCGIter,CGTol,x,xFeas);
        t = min([gamma, 1/(norm(d,'inf')^2)]);
        v = v + t*d;
        count = count + 1;
        CGIters(count) = cgIters_i;
        CGres(count)  = norm(r_i,2);
        resLim(count) = resLim_i;
        muVec(count) = mu;
        dNorm = norm(d,'inf');
        if printFlag
            fprintf('mu = %0.2e, d = %0.4f \n',mu,norm(d,'inf'))
        end
    end
end

% Check primal-dual feasibility condition
p = (Aineq')\r_i;
lambda = sqrt(mu)*(exp(v) + exp(v).*d) + p;
if min(lambda) > -1e-6
    feasFlag = 1;
else
    feasFlag = 0;

    % Run the post-processing step
    [x,d,cg] = solveNewtonStep(mu,v,const,maxCGIter,CGTol,x,x,1);
    alpha = min(1, 1/(norm(d,'inf')^2));
    v = v + alpha*d;
    CGIters(count) = CGIters(count) + cg;
    if printFlag
        fprintf('mu = %0.2e, d = %0.4f (Post-process) \n',mu,norm(d,'inf'))
    end



    while norm(d,'inf') > 1
        count = count + 1;
        % Run the post-processing step
        [x,d,cg] = solveNewtonStep(mu,v,const,maxCGIter,CGTol,x,x,1);
        alpha = min(1, 1/(norm(d,'inf')^2));
        v = v + alpha*d;
        CGIters(count) = cg;
        if printFlag
            fprintf('mu = %0.2e, d = %0.4f (Post-process) \n',mu,norm(d,'inf'))
        end
    end
end


% Set output
output.v = v;
output.CGIters = CGIters(1:count,:);
output.CGres = CGres(1:count,:);
output.muVec = muVec(1:count);
output.resLim =  resLim(1:count);
output.numIter = count;
output.maxIter = maxIter;
output.feasFlag = feasFlag;
output.lambda = lambda;
end



function [x,d,numIter,r,resLim,xFeas] = solveNewtonStep(mu,v,const,maxIter,tol,x0,xFeas,exactFlag)
if nargin < 8
    exactFlag = 0;
end

W = const.W;
% invW = const.invW;
c = const.c;
A = const.A;
b = const.b;
m = size(A,1);
gamma = const.gamma;
minEig = const.minEig;

% Define the RHS vector b
Q = diag(exp(2*v));
M = A'*Q*A + W;
expv = exp(v);
f = 2*sqrt(mu)*A'*expv - (c + A'*(exp(2*v).*b));

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

% Compute slack and Newton step
s_i = A*x+b;
d = ones(m,1) - 1/sqrt(mu)*expv.*s_i;

% Check if x is primal-feasible, otherwise use xFeas
if min(s_i) > 0
    xFeas = x;
    feasFlag = 1;
else
    feasFlag = 0;
end

% Then, depending on feasFlag, compute the upper-bound
if feasFlag
    sInv_i = 1./s_i;
    barGrad_i = A'*sInv_i;
    xDiff_1 = (2/minEig)*norm(W*x + c -  mu*barGrad_i,2);
    sigma = xDiff_1/mu;
else
    sInv_i = 1./(A*xFeas + b);
    barGrad_i = A'*sInv_i;
    xDiff_1 = (2/minEig)*norm(W*x + c -  mu*barGrad_i,2);
    xDiff_2 = norm(x - xFeas,2);
    sigma = (1/mu)*(xDiff_1 + xDiff_2);
end
% Check truncation criteria
dNorm = norm(d,2);
cond1 = (1-gamma)*dNorm^2;
resLim = cond1/sigma;
if res < resLim && ~exactFlag
    truncCond = 1;
else
    truncCond = 0;
end

% Iterate...
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

    % Compute slack and Newton step
    s_i = A*x+b;
    d = ones(m,1) - 1/sqrt(mu)*expv.*s_i;

    % Check if x is primal-feasible, otherwise use xFeas
    if min(s_i) > 0
        xFeas = x;
        feasFlag = 1;
    else
        feasFlag = 0;
    end

    % Then, depending on feasFlag, compute the upper-bound
    if feasFlag
        sInv_i = 1./s_i;
        barGrad_i = A'*sInv_i;
        xDiff_1 = (2/minEig)*norm(W*x + c -  mu*barGrad_i,2);
        sigma = xDiff_1/mu;
    else
        % Use the previous value of xDiff_1
        xDiff_2 = norm(x - xFeas,2);
        sigma = (1/mu)*(xDiff_1 + xDiff_2);
    end

    % Check truncation criteria
    dNorm = norm(d,2);
    cond1 = (1-gamma)*dNorm^2;
    resLim = cond1/sigma;
    if res < resLim && ~exactFlag
        truncCond = 1;
    end
end
end

function z = MTimes(x,v,A,W)

end
