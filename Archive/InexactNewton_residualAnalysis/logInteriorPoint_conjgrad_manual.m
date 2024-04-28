function [x,output] = logInteriorPoint_conjgrad_analysis(W,c,A,b,v0,opts)
% min 0.5*x'*W*x + c'*x   subject to:  A*x <= b

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
if isfield(opts,'maxIter')
    maxIter = opts.maxIter;
else
    maxIter = 150;
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

% Get size variables
m = size(A,1);
n = size(A,2);

% Store CG output feedback
x = zeros(n,1);

% First, change variables to Ax + b >= 0... This is just for uniformity
% with quadprog's inputs.
A = -A;

% Pack
invW = inv(W);
const.W = W;
const.invW = invW;
const.c = c;
const.A = A;
const.b = b;

% Initialize storage variables
pertDist_v = zeros(maxIter,1);
pertDist_x = zeros(maxIter,1);
distToCent_v = zeros(maxIter,1);
distToCent_x = zeros(maxIter,1);
distToPert_v = zeros(maxIter,1);
distToPert_x = zeros(maxIter,1);
dNorm = zeros(maxIter,1);
dTrueNorm = zeros(maxIter,1);
dError = zeros(maxIter,1);
resNorm = zeros(maxIter,1);
CGIters = zeros(maxIter,1);
% --------------------- NEWTON ITERATION ---------------------
% Initialize
v = v0;
mu = mu_0;
ACon = const.A;
bCon = const.b;

% Unpack centered point
x_c = opts.x_c;
v_c = opts.v_c;

% Iterate...
const_i = const;
for i = 1:maxIter
    % Define the RHS vector b
    f = ones(m,1) - 1/sqrt(mu)*exp(v).*(ACon*invWTimes(sqrt(mu)*ACon'*exp(v) - c,const) + bCon);

    % Get true Newton step for comparison
    [x_true,d_true] = solveNewtonStep(mu,v,const);

    % Perturb the Newton step by the true value
    randVec =  -1 + 2*rand(length(b),1); % random vector between -1 and 1
    r = (CGTol/norm(randVec))*randVec;

    % Get "d" for the perturbed point
    bTilde = b + r;
    const_i.b = bTilde;
    [x,d] = solveNewtonStep(mu,v,const_i);

    % Get the perturbed centered point
    centTol = 1e-6;
    [x_r,v_r,d_r] = logInteriorPoint_getCenteredPoint(W,c,-A,bTilde,v,mu,150,centTol);
    if norm(d_r,'inf') < centTol
        centFlag = 1;
    else
        centFlag = 0;
    end
    
    % Undo the change of variables and update v
    alpha_step = min(1, 1/(norm(d,'inf')^2));
    v = v + alpha_step*d;

    % Get distances and store
    if centFlag
        pertDist_v(i) = norm(v_c-v_r);
        pertDist_x(i) = norm(x_c-x_r);
        distToPert_v(i) = norm(v - v_r);
        distToPert_x(i) = norm(x - x_r);
    else
        pertDist_v(i) = nan;
        pertDist_x(i) = nan;
        distToPert_v(i) = nan;
        distToPert_x(i) = nan;
    end
    distToCent_x(i) = norm(x - x_c);
    distToCent_v(i) = norm(v - v_c);
    dNorm(i) = norm(d,2);
    dTrueNorm(i) = norm(d_true,2);
    dError(i) = norm(d - d_true);
    resNorm(i) = norm(r);
    CGIters(i) = 1;

    % Break condition
    fprintf('Error: %0.2e \n',distToCent_x(i))
    if distToCent_x(i) < mu_f
        break
    end
end

% Solve for primal variable x
output.v = v;

% Set outputs
maxIter = i;
output.maxIter = maxIter;
output.CGIters = CGIters(1:maxIter);
output.vr_to_vc = pertDist_v(1:maxIter);
output.xr_to_xc = pertDist_x(1:maxIter);
output.v_to_vc = distToCent_v(1:maxIter);
output.x_to_xc = distToCent_x(1:maxIter);
output.v_to_vr = distToPert_v(1:maxIter);
output.x_to_xr = distToPert_x(1:maxIter);
output.dNorm = dNorm(1:maxIter); 
output.dTrueNorm = dTrueNorm(1:maxIter);
output.dError = dError(1:maxIter);
output.resNorm = resNorm(1:maxIter); 
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


function [x,d] = solveNewtonStep(mu,v,const)
W = const.W;
c = const.c;
A = const.A;
b = const.b;
m = size(A,1);
n = size(A,2);

% Construct Q matrix
Qvec = exp(2*v);
Q = diag(Qvec);

% Solve for x
x = (A'*Q*A + W)\(2*sqrt(mu)*A'*exp(v) - (c + A'*Q*b));

% Solve for d
d = ones(m,1) - (1/sqrt(mu))*exp(v).*(A*x + b);

end


