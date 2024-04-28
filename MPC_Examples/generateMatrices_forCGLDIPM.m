function [H,F,M,L,b,AHat,BHat] = generateMatrices_forCGLDIPM(N,A,B,P,Q,R,xBar,uBar)
% * Generates the matrices H, g, Ain, bin to use in a quadratic program
% * from linear state-space MPC parameters. Unconstrained output and
% * control elements can be specified as 'NaN', these do NOT need to be
% * identical between ymax and ymin or umax and umin.

% Size variables
n = size(A,1);
m = size(B,2);

% First, construct the intermediate cost function matrices AA, BB, QQ, RR
AHat = zeros(n*N,n);
for i = 1:N
    AHat(1+(i-1)*n : n+(i-1)*n,:) = A^i;
end

BHat = zeros(n*N,m*N);
for i = 1:N % Skip the first row since its zeros
    for j = 1:i
        BHat(1+(i-1)*n : n+(i-1)*n, 1+(j-1)*m : m+(j-1)*m) = A^(i-j)*B;
    end
end

% COST MATRICES 

% intermediate matrices 
HHat = blkdiag(kron(eye(N-1),Q),P);
RHat = kron(eye(N),R);

% Now, construct the final Hessian matrix H and vector g
H = (RHat + BHat'*HHat*BHat);
[H] = conditionHessianMatrix(H);
F = BHat'*HHat*AHat;

% Next, generate the associated  constraint matrices
% st the constraints are represted by ACon*U <= FCon + LCon*xk

% Intermediate CHat and DHat matrices
EHat_x  = kron(eye(N),[eye(n); -eye(n)]);
EHat_u  = kron(eye(N),[eye(m); -eye(m)]);
CHat = [EHat_x; zeros(size(EHat_u,1),size(EHat_x,2))];
DHat = [zeros(size(EHat_x,1),size(EHat_u,2)); EHat_u];
b = [kron(ones(2*N,1),xBar);...
     kron(ones(2*N,1),uBar)];


%  Constraint matrice
M = CHat*BHat + DHat;
L = CHat*AHat;
