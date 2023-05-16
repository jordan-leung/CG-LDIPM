function [x,output] = logInteriorPoint_conjgrad_shortstep(W,c,A,b,v0,opts)
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
if isfield(opts,'beta')
    beta = opts.beta;  % size of the operational region 
else
    beta = 0.45; 
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
h_vHat_v = zeros(maxIter,1);
CGIters = zeros(maxIter,1);
muVec = zeros(maxIter,1);
resNorm = zeros(maxIter,1);

% Initialize
v = v0;
mu = mu_0;
ACon = const.A;
bCon = const.b;
centTol = 1e-6;

% Define shortstep parameters
epsilon = 0.5*qInv(beta);
a = beta; % init
count = 0;
while beta^(2*count) > epsilon^2
    count = count + 1;
end
N = count  % number of inner-loop iterations
zeta  = qInv(beta) - epsilon;
kappa = exp(2*qInv(zeta^2/m)) % 1/kappa reduction factor for mu

% ---- Begin by finding the centered point for the initial guess -------
% Get the centered point
[x_c,v_c] = logInteriorPoint_getCenteredPoint(W,c,-A,b,v,mu,150,centTol);
v = v_c;

% Main outer-loop for incrementing mu
iOuter = 1;
while mu > mu_f
    % Increment mu
    mu  = mu/kappa;

    % Get the centered point
    [x_c,v_c] = logInteriorPoint_getCenteredPoint(W,c,-A,b,v,mu,150,centTol);

    for i = 1:N % iterates for static mu
        % Define the RHS vector b
        f = ones(m,1) - 1/sqrt(mu)*exp(v).*(ACon*invWTimes(sqrt(mu)*ACon'*exp(v) - c,const) + bCon);

        % Initialize and redefine the problem such that x0 = 0.
        d0 = zeros(size(b,1),1); % correct at the end by d = x + d0
        d = d0;
        Md0 = MTimes(d,v,const);
        r = f - Md0;

        % Calculate iteration constants
        z = r;
        p = z;
        w = MTimes(p,v,const);
        alpha = r'*z/(p'*w);

        % Update x and r
        d = d + alpha*p;
        rPrev = r;
        r = r - alpha*w; % r1
        numIter_cg = 1;

        % --------------- CONJUGATE GRADIENT ---------------
        while numIter_cg  < maxCGIter && norm(r) > CGTol
            % Perform CG iterations
            zPrev = z;
            z = r;
            beta = r'*z/(rPrev'*zPrev);
            p = z + beta*p;
            w = MTimes(p,v,const);
            alpha = r'*z/(p'*w);

            % Update x
            d = d + alpha*p;
            rPrev = r;
            r = r - alpha*w;

            % i++
            numIter_cg = numIter_cg + 1;
        end

        % Undo the change of variables and update v
        d = d + d0;
        x = invWTimes(sqrt(mu)*A'*(exp(v) + exp(v).*d) - c,const);
        alpha_step = min(1, 1/(norm(d,'inf')^2));
        v = v + alpha_step*d;

        % Get distances and check the break criteria
        h_mu = hDiv(v_c,v);


        % Get distances and store
        h_vHat_v(iOuter) = h_mu;
        CGIters(iOuter) = numIter_cg;
        muVec(iOuter) =  mu;
        resNorm(iOuter) = norm(r);
%         fprintf('--- iter  = %0.0i, mu = %0.2e ---- \n',iOuter,mu)
        iOuter = iOuter + 1;
    end
end % end while loop

% Solve for primal variable x
output.v = v;

% Set outputs
maxIter = iOuter-1;
output.maxIter = maxIter;
output.CGIters = CGIters(1:maxIter);
output.h_vHat_v = h_vHat_v(1:maxIter);
output.muVec = muVec(1:maxIter);
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

function [h] = hDiv(v1,v2)
h = exp(v1)'*exp(-v2) + exp(-v1)'*exp(v2) - 2*length(v1);
end

function p = pFunc(t)
p = t - qFunc(qFunc(sqrt(t)));
end

function q = qFunc(t)
q = 2*(cosh(t) - 1);
end

function q = qInv(theta)
q = acosh(0.5*theta + 1);
end