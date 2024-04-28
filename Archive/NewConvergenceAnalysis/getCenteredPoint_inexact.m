function [x,v,output] = getCenteredPoint_base(W,c,A,b,x0,v0,mu,maxIter,dTol,v_c,gamma)
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
const.mm = min(eig(W));
const.gamma = gamma;
% Outputs
numIterVec = zeros(maxIter,1);
hVec = zeros(maxIter,1);
dNormVec = zeros(maxIter,1);
rNormVec = zeros(maxIter,1);
sigma = zeros(maxIter,1);
delta = zeros(maxIter,1);

% Get centered point
x = x0;
for i = 1:maxIter
    [x,d,numIter_i,rNorm_i,sigma_i,delta_i] = solveNewtonStep(mu,v,const,1000,1e-8,x);
    normd = norm(d,'inf');
    alpha = min(gamma, 1/(normd^2));
%         alpha = min(gamma, 1/(normd^2));
    vPrev = v;
    v = v + alpha*d;
    if normd < dTol
        break
    end
    numIterVec(i) = numIter_i;
    hVec(i) = hDiv(vPrev,v_c);
    dNormVec(i) = norm(d,2);
    rNormVec(i) = rNorm_i;
    sigma(i) = sigma_i;
    delta(i) = delta_i;
end
output.numIterVec = numIterVec(1:i);
output.hVec = hVec(1:i);
output.dNormVec = dNormVec(1:i);
output.rNormVec = rNormVec(1:i);
output.sigma = sigma(1:i);
output.delta = delta(1:i);
end


function [x,d,numIter,res,sigma,delta] = solveNewtonStep(mu,v,const,maxIter,tol,x0)
W = const.W;
% invW = const.invW;
c = const.c;
A = const.A;
b = const.b;
mm = const.mm;
m = size(A,1);
gamma = const.gamma;

% Define the RHS vector b
Q = diag(exp(2*v));
M = A'*Q*A + W;
f = 2*sqrt(mu)*A'*exp(v) - (c + A'*Q*b);

% --------------- CONJUGATE GRADIENT ---------------
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

% Check truncation criteria
sigma = (2/(mm*mu))*norm(W*x + c + mu*barrierGrad(x,A,b),2);
d = ones(m,1) - 1/sqrt(mu)*exp(v).*(A*x + b);
delta = (1-gamma)*norm(d,2)^2;
if sigma*res < delta
    truncCond = 1;
else
    truncCond = 0;
end
while numIter  < maxIter && res > tol && truncCond == 0
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
    sigma = (2/(mm*mu))*norm(W*x + c + mu*barrierGrad(x,A,b),2);
    d = ones(m,1) - 1/sqrt(mu)*exp(v).*(A*x + b);
    delta = (1-gamma)*norm(d,2)^2;
    if sigma*res < delta
        truncCond = 1;
    end
end
end

function [h] = hDiv(v1,v2)
h = exp(v1)'*exp(-v2) + exp(-v1)'*exp(v2) - 2*length(v1);
end

function val = barrierGrad(x,A,b)
   m = size(A,1);
   n = size(A,2);
   val = zeros(n,1);
   for i = 1:m
       a_i = A(i,:)';
       b_i = b(i);
       val = val - a_i/(a_i'*x+b_i);
   end
end
