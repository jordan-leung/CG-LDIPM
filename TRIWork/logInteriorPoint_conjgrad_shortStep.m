function [x,v,mu,execTime,numIter,CGIters,CGres,CGerror,dHist] = logInteriorPoint_conjgrad_shortStep(W,c,Aineq,bineq,mu_f,mu_0,v0,maxIter,maxCGIter,CGTol,preCondFlag,params,vThresh,vNumThresh)
% min 0.5*x'*W*x + c'*x   subject to:  A*x <= b
% Get size variabl,es
m = size(Aineq,1);
n = size(Aineq,2);

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
A = Aineq;

% Compute the diagonal values of G which we use for preconditioning in CG
GDiag = zeros(m,1);
const.GDiag = GDiag;
for i = 1:m
    a_i = A(i,:);
    GDiag(i) = a_i*invWTimes(a_i',const);
end
const.GDiag = GDiag;

% Max amount of bisection iterations
numIter = 0; % number of newton iterations performed

% Initialize things
CGIters = zeros(maxIter,1);
CGres = zeros(maxIter,1);
CGerror = zeros(maxIter,1);
dHist = zeros(maxIter,1);

% Given parameters (delta,epsilon), deteremine the shortstep parameters (N,k)
N_ls = params.N;
k_ls = params.k;
const.params = params;
const.vThresh = vThresh;
const.vNumThresh = vNumThresh;

tic
% --------------------- INITAL CENTERING PROCEDURE ---------------------
dNorm = 1;
v = v0;
mu = mu_0;
d = zeros(m,1);
while dNorm > 1e-8
    [d,CGIter_i,res] = solveNewtonStep(mu,v,const,d,preCondFlag,maxCGIter,CGTol);
    dNorm = norm(d,'inf');
    fprintf('mu = %0.2e, d = %0.4f (Centering) \n',mu,dNorm)
    
    % For calculating error
    MTemp = eye(m) + diag(exp(v))*A*invW*A'*diag(exp(v));
    f = ones(m,1) - 1/sqrt(mu)*exp(v).*(A*invWTimes(sqrt(mu)*A'*exp(v) - c,const) + bineq);
    dOpt = MTemp\f;
    e_i = norm(d - dOpt,2);
        
    
    % Update x, v, d
    alpha = min(1, 1/(dNorm^2));
    v = v + alpha*d;
    numIter = numIter + 1;
    
    % Store
    CGIters(numIter) = CGIter_i;
    CGres(numIter) = res;
    CGerror(numIter) = e_i;
    dHist(numIter) = dNorm;
end
% --------------------- MAIN NEWTON ITERATION LOOP ---------------------
dPrev = d;
while mu > mu_f
    % Perturb mu
    mu = (1/k_ls)*mu;

    % Run N inner-loop iterations
    for j = 1:N_ls
        % Run the Newton system.
        [d,CGIter_i,res] = solveNewtonStep(mu,v,const,d,preCondFlag,maxCGIter,CGTol);
        
        % For calculating error
        MTemp = eye(m) + diag(exp(v))*A*invW*A'*diag(exp(v));
        f = ones(m,1) - 1/sqrt(mu)*exp(v).*(A*invWTimes(sqrt(mu)*A'*exp(v) - c,const) + bineq);
        dOpt = MTemp\f;
        e_i = norm(d - dOpt,2);
        
        % Update x, v, d
        dNorm = norm(d,'inf');
        vPrev = v;
        alpha = min(1, 1/(dNorm^2));
        v = v + alpha*d;
        fprintf('mu = %0.2e, d = %0.4f \n',mu,dNorm)
        
        % Store
        numIter = numIter + 1;
        CGIters(numIter) = CGIter_i;
        CGres(numIter) = res;
        CGerror(numIter) = e_i;
        dHist(numIter) = norm(d,2);
    end
end
execTime = toc;
v = vPrev;
x = invW*(sqrt(mu)*A'*(exp(v) + exp(v).*d) - c);
CGIters = CGIters(1:numIter);
CGres = CGres(1:numIter);
CGerror = CGerror(1:numIter);
dHist = dHist(1:numIter);
end


% This is a placeholder function for when we eventually use Riccatti
function zOut = invWTimes(zIn,const)
% Returns zOut = invW*zIn
zOut = const.invW*zIn;
end

% Function to evaluate M(v)*x
function zOut = MTimes(zIn,v,const)
% Returns zOut = M(v)*zIn
A = const.A;
D = exp(v);
zOut = zIn + D.*(A*( invWTimes(A'*(D.*zIn),const)));
end


function [d,numIter,res] = solveNewtonStep(mu,v,const,d0,preCondFlag,maxIter,tol)
% W = const.W;
% invW = const.invW;
c = const.c;
ACon = const.A;
bCon = const.b;
m = size(ACon,1);
GDiag = const.GDiag;

% Define the preconditioner MTilde
MTilde = ones(m,1) + exp(v).*GDiag.*exp(v);

% Define the RHS vector b
f = ones(m,1) - 1/sqrt(mu)*exp(v).*(ACon*invWTimes(sqrt(mu)*ACon'*exp(v) - c,const) + bCon);

% % % For debugging purposes 
% MTemp = eye(m) + diag(exp(v))*ACon*const.invW*ACon'*(exp(v));
% dOpt = MTemp\f;

% First, determine whether or not to apply the diagonal preconditioner. Use
% a criterion than at least 1/4 of the variables have dropped below
% vTresh... This means that many elements of exp(v) will be near zero
vThresh = const.vThresh;
vNumThresh = const.vNumThresh;
negVec = find(v < vThresh);
if preCondFlag == 1
    if length(negVec) > vNumThresh
        applyPreCond = 1;
    else
        applyPreCond = 0;
    end
else
    applyPreCond = 0;
end

% --------------- CONJUGATE GRADIENT ---------------
% Initialize and redefine the problem such that x0 = 0.
x = d0; % correct at the end by d = x + d0
Md0 = MTimes(d0,v,const);
b = f - Md0;

% Run the first iteration of CG and iniialize iteration variables
r = b;

% Calculate iteration constants
if applyPreCond
    z = r./MTilde; % preconditioner step
else
    z = r;
end
p = z;
w = MTimes(p,v,const);
alpha = r'*z/(p'*w);

% Update x and r
x = x + alpha*p; 
rPrev = r;
r = r - alpha*w; % r1
res = norm(r,2);
numIter = 1;

% Iterate... 
dNorm = norm(x);
while numIter  <= maxIter && res > tol && ~(checkTruncCriteria(res,dNorm,const.params))
    zPrev = z;
    if applyPreCond
        z = r./MTilde;
    else
        z = r;
    end
    beta = r'*z/(rPrev'*zPrev);
    p = z + beta*p;
    w = MTimes(p,v,const);
    alpha = r'*z/(p'*w);
    
    % Update x
    x = x + alpha*p;
    rPrev = r;
    r = r - alpha*w;
    res = norm(r,2);

    % Update the dNorm
    dNorm = norm(x);
    
    % i++
    numIter = numIter + 1;
end
d = x;

end

function flag = checkTruncCriteria(eNorm,dNorm,params)
% Returns 1 if the CG truncation criteria is satisifed 
delta = params.delta;
gamma = delta/2 + sqrt(delta^2/4 + delta);

% Compute upper and lower bounds of h
d_ub = min([dNorm + eNorm, gamma]);
h_ub = min([(d_ub)^2/(1 - dNorm - eNorm), delta]);
h_lb = (dNorm - eNorm)^2/(1 + d_ub);

% Compute the complicated LHS and RHS for the 3 bigger inequalities
LHS_delta = delta/2 + eNorm + sqrt(delta^2/4 + 2*delta*eNorm + delta);
rho = (2*dNorm + sqrt(h_ub^2 + 4*h_ub))*eNorm + eNorm^2;
LHS_quad = rho + h_ub*dNorm*eNorm;
RHS_quad = 0.5*(1 - dNorm - eNorm)*h_lb^2;

if eNorm > dNorm || eNorm > 1 - dNorm || LHS_delta > 1 || LHS_quad > RHS_quad
    flag = 0;
else
    flag = 1;
end
end