function [x,v,output] = getCenteredPoint_base(W,c,A,b,x0,v0,mu,maxIter,dTol,v_c)
% min 0.5*x'*W*x + c'*x   subject to:  A*x <= b


% Get size variables
m = size(A,1);
n = size(A,2);

% First, change variables to Ax + b >= 0... This is just for uniformity
% with quadprog's inputs.
A = -A;
v = v0;

% Pack
invW = inv(W);
const.W = W;
const.invW = invW;
const.c = c;
const.A = A;
const.b = b;

% Outputs
numIterVec = zeros(maxIter,1);
hVec = zeros(maxIter,1);
dNormVec = zeros(maxIter,1);

% Get centered point
x = zeros(n,1);
for i = 1:maxIter
    [x,d,numIter_i] = solveNewtonStep(mu,v,const,1000,1e-8,zeros(n,1));
    normd = norm(d,'inf');
    alpha = min(1, 1/(normd^2));
    vPrev = v;
    v = v + alpha*d;
    if normd < dTol
        break
    end
    numIterVec(i) = numIter_i;
    hVec(i) = hDiv(vPrev,v_c);
    dNormVec(i) = norm(d,2);
end
output.numIterVec = numIterVec(1:i);
output.hVec = hVec(1:i);
output.dNormVec = dNormVec(1:i);
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
d = 1 - 1/sqrt(mu)*exp(v).*(A*x + b);
end


function [h] = hDiv(v1,v2)
h = exp(v1)'*exp(-v2) + exp(-v1)'*exp(v2) - 2*length(v1);
end
