function [fHist,dHist,fPrimeHist,fPrimeInitHist,CGIters,CGIters_thresh,CGIters_decrease] = divergenceAnalysis(W,c,Aineq,bineq,mu,v0,CGTol,preCondFlag,vThresh,vNumThresh,numSamples)
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
const.vThresh = vThresh;
const.vNumThresh = vNumThresh;


% Initialize things
dHist = zeros(numSamples,1);
fHist = zeros(numSamples,1);
fPrimeHist = zeros(numSamples,1);
fPrimeInitHist = zeros(numSamples,1);
CGIters = zeros(numSamples,1);
CGIters_thresh = zeros(numSamples,1);
CGIters_decrease = zeros(numSamples,1);

% --------------------- INITAL CENTERING PROCEDURE ---------------------
dNorm = 1;
v = v0;
d = zeros(m,1);
while dNorm > 1e-10
    [d,CGIter_i,res] = solveNewtonStep(mu,v,const,d*0,preCondFlag,1e8,CGTol,[]);
    dNorm = norm(d,'inf');
    fprintf('mu = %0.2e, d = %0.4f (Centering) \n',mu,dNorm)
    
    % Update x, v, d
    alpha = min(1, 1/(dNorm^2));
    v = v + alpha*d;
end


% --------------------- MAIN SAMPLING LOOP ---------------------
vHat = v;
sig_max = 5*norm(vHat,2);
sig_min = 1e-4*norm(vHat,2);
for i = 1:numSamples
    % Generate random vector to perturb v
    sig_i = sig_min * (sig_max-sig_min)*rand(1,1);
    v = vHat + normrnd(zeros(m,1),sig_i*ones(m,1));

    % Calculate Newton vector
    [d,CGIter_i,res,theshIter_i,decreaseIter_i,fPrimeInit_i] = solveNewtonStep(mu,v,const,d*0,preCondFlag,1e8,CGTol,vHat);
    
    % Calculate divergence
    h_i = (exp(vHat))'*(exp(-v)) + (exp(-vHat))'*(exp(v)) - 2*m;
    fPrime_i = (exp(v - vHat) - exp(vHat - v))'*d;
    dNorm = norm(d,2);
    
    % Store
    fHist(i) = h_i;
    dHist(i) = dNorm;
    fPrimeHist(i) = fPrime_i;
    fPrimeInitHist(i) = fPrimeInit_i;    
    CGIters(i) = CGIter_i;
    CGIters_thresh(i) = theshIter_i;
    CGIters_decrease(i) = decreaseIter_i;
end
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


function [d,numIter,res,numIter_thresh,numIter_decrease,fPrimeInit] = solveNewtonStep(mu,v,const,d0,preCondFlag,maxIter,tol,vHat)
% W = const.W;
% invW = const.invW;
c = const.c;
ACon = const.A;
bCon = const.b;
m = size(ACon,1);
GDiag = const.GDiag;

% Flag for whether or not we use the divergence criteria
checkDiv = ~isempty(vHat);
numIter_thresh = 0;
numIter_decrease = 0; % assign for cases where we don't use this
fPrimeInit = 0;

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

% Calculate initial fPrime
if checkDiv
    fPrimeInit = (exp(v - vHat) - exp(vHat - v))'*x;
end

% Iterate... 
truncFlag_flag = 0;
decreaseFlag_flag = 0;
% Check truncation criterion
if checkDiv
    [truncFlag,decreaseFlag] = checkTruncCriteria(x,v,vHat); % note that x = d (candidate) here
    
    % Check if the "thresh" condition has been met
    if ~truncFlag_flag
        if truncFlag == 1
            numIter_thresh = numIter;
            truncFlag_flag = 1;
        end
    end
    
    % Check if f'(0) has been met
    if ~decreaseFlag_flag
        if decreaseFlag == 1
            numIter_decrease = numIter;
            decreaseFlag_flag = 1;
        end
    end
end


while numIter  <= maxIter && res > tol % go until the tolerance is hit, but check for indices
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
        [truncFlag,decreaseFlag] = checkTruncCriteria(x,v,vHat); % note that x = d (candidate) here
        
        % Check if the "thresh" condition has been met
        if ~truncFlag_flag
            if truncFlag == 1
                numIter_thresh = numIter;
                truncFlag_flag = 1;
            end
        end
        
        % Check if f'(0) has been met
        if ~decreaseFlag_flag
            if decreaseFlag == 1
                numIter_decrease = numIter;
                decreaseFlag_flag = 1;
            end
        end
    end
end
d = x;
end

function [flag,flag2] = checkTruncCriteria(d,v,vHat)
% Returns 1 if the CG truncation criteria is satisifed 
f = (exp(vHat))'*(exp(-v)) + (exp(-vHat))'*(exp(v)) - 2*length(v);
fPrime = (exp(v - vHat) - exp(vHat - v))'*d;
if fPrime < -f - norm(d,2)^2
    flag = 1;
else
    flag = 0;
end

if fPrime < 0
    flag2 = 1;
else
    flag2 = 0;
end

end