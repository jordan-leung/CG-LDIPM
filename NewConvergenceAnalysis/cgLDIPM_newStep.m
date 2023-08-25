function [x,output] = cgLDIPM_newStep(W,c,Aineq,bineq,v0,opts)
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
% invW = inv(W);
const.W = W;
% const.invW = invW;
const.c = c;
const.A = Aineq;
const.b = bineq;
const.minEig = min(eig(const.W));
numIter = 0; % number of newton iterations performed

% --------------------- MAIN ITERATION LOOP ---------------------
% Unpack shortstep parameters
v = v0;
x = zeros(n,1);
mu = mu_0;
kappa  =  opts.kappa;
N = opts.N;
gamma = opts.gamma;
const.gamma = gamma;
const.minEig = min(eig(W));

% Store CG output feedback
maxIter = ceil(1*log(mu_0/mu_f)*sqrt(m)/(2*qInv(opts.zeta^2))); % just for initialization
CGIters = zeros(maxIter,1);
CGres = zeros(maxIter,1);
muVec = zeros(maxIter,1);
count = 0;
while mu > mu_f
    % Compute d
    [x,d,cgIters_i,res_i,] =  solveNewtonStep(mu,v,const,maxIter,CGTol,x);

    % Update mu and d
    if count > 500
        count
    end
    [d,mu] = muUpdate(mu,x,v,d,const);
    
    % Step
    t = min([gamma, 1/(norm(d,'inf')^2)]);
    v = v + t*d;

    % Increment
    count = count + 1;
    CGIters(count) = cgIters_i;
    CGres(count)  = res_i;
    muVec(count) = mu;
    fprintf('mu = %0.2e, d = %0.4f  \n',mu,norm(d,'inf'))
end

% Set output
output.v = v;
output.CGIters = CGIters(1:count,:);
output.CGres = CGres(1:count,:);
output.muVec = muVec(1:count);
output.numIter = count;
output.maxIter = maxIter;
end



function [x,d,numIter,res] = solveNewtonStep(mu,v,const,maxIter,tol,x0)
W = const.W;
% invW = const.invW;
c = const.c;
A = const.A;
b = const.b;
m = size(A,1);

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
while numIter  < maxIter && res > tol 
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
end
d = ones(m,1) - 1/sqrt(mu)*exp(v).*(A*x + b);
end




function [d,mu] = muUpdate(mu0,x,v,d_in,const)


% Unpack
A = const.A;
b = const.b;
W = const.W;
c = const.c;
gamma = const.gamma;
m = size(A,1);

% Initial candidate mu
mu_true = mu0;
d_true = d_in;
maxIter = 1000; % maximum bisection iterations

% Pre-calculate
expv = exp(v);
s = A*x + b;
sInv = 1./s;
barGrad = A'*sInv;
Wx = W*x;
term1 =  A'*((expv.^2).*s);
oneVec = ones(m,1);

% Execute decreases until we find a non-valid mu
numIter = 0;
mu = mu0;
while true
    % Incrmeent mu
    numIter = numIter + 1;
    mu = 0.9*mu;

    % Update vectors
    dPlus = 1 - 1/sqrt(mu)*expv.*s;
    rPlus = (sqrt(mu0/mu) - 1)*term1 +  (sqrt(mu0) - sqrt(mu)) + A'*(expv + expv.*dPlus);
    
    % Check truncation criteria
    sigma = (2/(const.minEig*mu))*norm(Wx + c -  mu*barGrad,2);
    if sigma*norm(rPlus) > (1-gamma)*norm(dPlus)
        mu_false = mu;
        break
    else
        mu_true = mu;
        d_true = dPlus;
    end
end

% Run bisection
while numIter < maxIter || (mu_true - mu_false) > 0.01*mu
    % Set mu
    numIter = numIter + 1;
    mu = 0.5*(mu_true + mu_false);

    % Update vectors
    dPlus = 1 - 1/sqrt(mu)*expv.*s;
    rPlus = (sqrt(mu0/mu) - 1)*term1 +  (sqrt(mu0) - sqrt(mu)) + A'*(expv + expv.*dPlus);

    % Check truncation criteria
    sigma = (2/(const.minEig*mu))*norm(Wx + c -  mu*barGrad,2);
    if sigma*norm(rPlus) > (1-gamma)*norm(dPlus)
        mu_false = mu;
    else
        mu_true = mu;
        d_true = dPlus;
    end
end
mu = mu_true;
d = d_true;
end
