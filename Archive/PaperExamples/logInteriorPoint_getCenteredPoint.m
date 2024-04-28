function [x,v,d,v_init] = logInteriorPoint_getCenteredPoint(W,c,A,b,v0,mu,maxIter,dTol)
% min 0.5*x'*W*x + c'*x   subject to:  A*x <= b


% Get size variables
m = size(A,1);
n = size(A,2);

% First, run regular LDIPM to get close to the centered point
opts.mu_0 = 1e8;
opts.mu_f = mu;
opts.printFlag = 0;
opts.maxIter = 100;
[~,output_init] = logInteriorPoint(W,c,A,b,v0,opts);

% First, change variables to Ax + b >= 0... This is just for uniformity
% with quadprog's inputs.
A = -A;
v_init = output_init.v;
v = v_init;

% Pack
invW = inv(W);
const.W = W;
const.invW = invW;
const.c = c;
const.A = A;
const.b = b;

% Get centered point
for i = 1:maxIter
    [x,d] = solveNewtonStep(mu,v,const);
    normd = norm(d,'inf');
    alpha = min(1, 1/(normd^2));
    v = v + alpha*d;
    if normd < dTol
        break
    end
end


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

