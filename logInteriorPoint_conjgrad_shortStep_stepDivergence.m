function [x,v,mu,execTime,numIter,CGIters,CGres,dHist,divHist,decreaseIters] = logInteriorPoint_conjgrad_calcDivergence(W,c,Aineq,bineq,mu_f,mu_0,v0,maxIter,maxCGIter,CGTol,preCondFlag,params,vThresh,vNumThresh)
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

numIter = 0; % number of newton iterations performed
N_ls = params.N;
k_ls = params.k;

% Initialize things
CGIters = zeros(maxIter,1);
CGres = zeros(maxIter,1);
dHist = zeros(maxIter,1);
divHist = zeros(maxIter,1);
decreaseIters = zeros(maxIter,1);
const.vThresh = vThresh;
const.vNumThresh = vNumThresh;

tic
% --------------------- INITAL CENTERING PROCEDURE ---------------------
dNorm = 1;
v = v0;
mu = mu_0;
d = zeros(m,1);
while dNorm > 1e-8
    [d,CGIter_i,res] = solveNewtonStep(mu,v,const,d*0,preCondFlag,maxCGIter,CGTol,[]);
    dNorm = norm(d,'inf');
    fprintf('mu = %0.2e, d = %0.4f (Centering) \n',mu,dNorm)
    
    % Update x, v, d
    alpha = min(1, 1/(dNorm^2));
    v = v + alpha*d;
    numIter = numIter + 1;
    
    % Store
    CGIters(numIter) = CGIter_i;
    CGres(numIter) = res;
    dHist(numIter) = dNorm;
    decreaseIters(numIter) = Inf;
end
% --------------------- MAIN NEWTON ITERATION LOOP ---------------------
while mu > mu_f
    % Perturb mu
    mu = (1/k_ls)*mu;

    % Run to get centered point so we can cheat and use the divergence
    % explicitly
    vHat = v;
    while dNorm > 1e-10
        % Run the Newton system.
        [d,CGIter_i,res] = solveNewtonStep(mu,vHat,const,d*0,preCondFlag,maxCGIter,CGTol,[]);
        
        % Update x, v, d
        dNorm = norm(d,'inf');
        alpha = min(1, 1/(dNorm^2));
        vHat = vHat + alpha*d;
    end
    
    
    % Run inner-loop iterations until the divergence condition is
    % satisfied. Run outer-loop iterations until ||d|| < dSize (e.g. 0.5)
    d = d*0;
    dNorm = 1;
    for j = 1:N_ls
        % Run the Newton system.
        [d,CGIter_i,res,decreaseIter_i] = solveNewtonStep(mu,v,const,d*0,preCondFlag,maxCGIter,CGTol,vHat);
        
        % Calculate divergence
        h_i = (exp(vHat))'*(exp(-v)) + (exp(-vHat))'*(exp(v)) - 2*m; 

        % Update x, v, d
        dNorm = norm(d,'inf');
        alpha = min(1, 1/(dNorm^2));
        v = v + alpha*d;
        fprintf('mu = %0.2e, d = %0.4f \n',mu,dNorm)
                
        % Store
        numIter = numIter + 1;
        CGIters(numIter) = CGIter_i;
        CGres(numIter) = res;
        dHist(numIter) = dNorm;
        decreaseIters(numIter) = decreaseIter_i;
        divHist(numIter) = h_i;
    end
end
execTime = toc;
x = invW*(sqrt(mu)*A'*(exp(v) + exp(v).*d) - c);
CGIters = CGIters(1:numIter);
CGres = CGres(1:numIter);
dHist = dHist(1:numIter);
decreaseIters = decreaseIters(1:numIter);
divHist = divHist(1:numIter);
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


function [d,numIter,res,fPrime] = solveNewtonStep(mu,v,const,d0,preCondFlag,maxIter,tol,vHat)
% W = const.W;
% invW = const.invW;
c = const.c;
ACon = const.A;
bCon = const.b;
m = size(ACon,1);
GDiag = const.GDiag;

% Flag for whether or not we use the divergence criteria
checkDiv = ~isempty(vHat);

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
truncFlag = 0;
% Check truncation criterion
if checkDiv
    [truncFlag,fPrime] = checkTruncCriteria(x,v,vHat); % note that x = d (candidate) here
end
while numIter  <= maxIter && res > tol && truncFlag == 0
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

    % i++
    numIter = numIter + 1;
    
    % Check truncation criterion
    if checkDiv
        [truncFlag,fPrime] = checkTruncCriteria(x,v,vHat); % note that x = d (candidate) here
    end
end
d = x;

end

function [flag,fPrime] = checkTruncCriteria(d,v,vHat)
% Returns 1 if the CG truncation criteria is satisifed 
f = (exp(vHat))'*(exp(-v)) + (exp(-vHat))'*(exp(v)) - 2*length(v);
fPrime = (exp(v - vHat) - exp(vHat - v))'*d;
if fPrime < -f - norm(d,2)^2
    flag = 1;
else
    flag = 0;
end


end