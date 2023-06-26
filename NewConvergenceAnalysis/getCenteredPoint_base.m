function [x,v,output] = getCenteredPoint_base(W,c,A,b,v0,mu,maxIter,dTol,v_c)
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
for i = 1:maxIter
    [d,numIter_i] = solveNewtonStep(mu,v,const,1000,1e-8,zeros(m,1));
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
x = invWTimes(sqrt(mu)*A'*(exp(vPrev) + exp(vPrev).*d) - c,const);
output.numIterVec = numIterVec(1:i);
output.hVec = hVec(1:i);
output.dNormVec = dNormVec(1:i);
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


function [d,numIter] = solveNewtonStep(mu,v,const,maxIter,tol,d0)
% W = const.W;
% invW = const.invW;
c = const.c;
ACon = const.A;
bCon = const.b;
m = size(ACon,1);

% Define the RHS vector b
f = ones(m,1) - 1/sqrt(mu)*exp(v).*(ACon*invWTimes(sqrt(mu)*ACon'*exp(v) - c,const) + bCon);

% --------------- CONJUGATE GRADIENT ---------------
% Initialize and redefine the problem such that x0 = 0.
x = zeros(size(d0,1),1); % correct at the end by d = x + d0
Md0 = MTimes(d0,v,const);
b = f - Md0;

% Run the first iteration of CG and iniialize iteration variables
r = b;
xStar = x; % initialize
resStar = norm(r,2); % initialize

% Calculate iteration constants
z = r;
p = z;
w = MTimes(p,v,const);
alpha = r'*z/(p'*w);

% Update x and r
x = x + alpha*p; 
rPrev = r;
r = r - alpha*w; % r1
res = norm(r,2);
if res < resStar % Store d as minimal residual solution
    xStar = x;
    resStar = res;
end
numIter = 1;

% Iterate... 
while numIter  < maxIter && res > tol
    zPrev = z;
    z = r;
    beta = r'*z/(rPrev'*zPrev);
    p = z + beta*p;
    w = MTimes(p,v,const);
    alpha = r'*z/(p'*w);
    
    % Update x
    x = x + alpha*p;
    rPrev = r;
    r = r - alpha*w;
    res = norm(r,2);
    if res < resStar
        xStar = x;
        resStar = res;
    end

    % i++
    numIter = numIter + 1;
end

% Undo the change of variables
d = xStar + d0;
end

function [h] = hDiv(v1,v2)
h = exp(v1)'*exp(-v2) + exp(-v1)'*exp(v2) - 2*length(v1);
end
