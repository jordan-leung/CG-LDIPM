clc
clear all
close all

clear
close all
clc

addpath('/Users/jordan/Documents/GitHub/optimizationFunctions');
addpath('/Users/jordan/Documents/GitHub/extraFunctions')

set(groot,'defaultFigureCreateFcn',@(fig,~)addToolbarExplorationButtons(fig))

%% * Set simulation parameters

% * Time vector
dt = 0.01;
t = 0:dt:3;
N = length(t);

MASS = 2041; % vehicle mass, kg
Izz = 4964; % yaw moment-of-intertia, kg/m^2
ell_f = 1.56; % front moment arm
ell_r = 1.64; % rear moment arm
C_alpha = 246994; % tire stiffness, N/rad
Vx = 30; % longitudinal velocity, m/s


% State Vector (X): s          = lateral position, (n)
%                   \psi       = lateral angle (rad)
%                   beta       = sideslip angle (rad), yDot/Ux
%                   omega      = yaw rate (rad/s), yDot/Ux

% Control Vector:   \delta     = steering angle (rad)

% Constrained output (y): alpha_f  = front slip angle, (rad)
%                         alpha_r  = rear slip angle, (rad)
%                         \delta   =  steering angle, (rad)

% Tracking output (z): s          = lateral position, (n)

% * Initial Conditions
X0 = [0; 0; 0; 0];
% X0 = [-4.9; 0; 0; 0];
n = length(X0);

% * Define commanded reference
nr = 1;



% Continuous-time dynamics
Ac = [0, Vx, Vx, 0;...
     0, 0, 0, 1;...
     0, 0, -2*C_alpha/(MASS*Vx), C_alpha*(ell_r - ell_f)/(MASS*Vx^2) - 1;...
     0, 0, C_alpha*(ell_r - ell_f)/Izz, -C_alpha*(ell_r^2 + ell_f^2)/(Izz*Vx)];
 
 Bc = [0; 0; C_alpha/(MASS*Vx); C_alpha*ell_f/Izz];

% Get discrete state-space matrices
m = size(Bc,2);
[A,B] = c2d(Ac,Bc,dt);


%% Generate tracking system matrices and constraint sets

% % constraints in the form Cx + Du \in Y... Relates states to slip angles
CPrime = [0, 0, -1, -ell_f/Vx;...
          0, 0, -1, ell_r/Vx;...
          0, 0, 0, 0];
DPrime = [1; 0; 1];


% Constraints. Ay <= b
ny = size(CPrime,1);
Ay = [eye(ny); -eye(ny)];
ymax = [8; 8; 30]*pi/180; % OUTPUT CONSTRAINTS
ymin = -ymax;
by = [ymax; -ymin];

% Tracking sttates in the form z = Ex + Fu
Es = [1 0 0 0];
F = 0;

[nx,nu] = size(B);

Z = [eye(nx)-A, B, zeros(nx,nu); Es, F, -eye(nr)]; 
null_Z = null(Z);

Gx = null_Z(1:nx,:);
Gu = null_Z(nx+1:nx+nu,:);
Gr = null_Z(nx+nu+1:end,:);

if size(Gr,1) == size(Gr,2)
    if det(Gr) ~= 0
        Gx = Gx*inv(Gr);
        Gu = Gu*inv(Gr);
        Gr = eye(nr);
    end
end



%% Simulation and MPC parameters
N_MPC = 66;  % 66 for standard mpc at x0 = 0
T_MPC = dt;

% * Set P, Q, R
Q = 0.1*eye(n);
Q(1,1) = 1;
% Q = Es'*Es; %  State error cost
R = 0.1; % Control effort cost
[K,P] = dlqr(A,B,Q,R);


if mod(T_MPC,dt) ~= 0
    error('Mismatch in divisibility of T and dt')
end

%++++++++++++++++++++++++ Generate Oinf +++++++++++++++++++++++++++++++++++
A_bar = A-B*K;
B_bar = B*K*Gx + B*Gu;
C_bar = CPrime-DPrime*K;
D_bar = DPrime*K*Gx + DPrime*Gu;
  

% Generate Oinf and flips column order from [r x] to [x r]
[T_inf, c_inf] = genOinf(A_bar, B_bar, C_bar, D_bar , Ay, by, 200,.01,1e-12); % generate Oinf
temp = T_inf(:,nr+1:end);
T_inf = [temp, T_inf(:,1:nr)]; 

% Generate QP Matrices
[H_MPC,W_MPC,M_MPC,L_MPC,b_MPC,AHat,BHat] = generateQPMatrices_tracking(N_MPC,A,B,CPrime,DPrime,P,Q,R,Ay,by,T_inf,c_inf,Gx,Gu);
H_QP = H_MPC;
A_QP = M_MPC;

       
%% Run Simulation

% Initialize trajectories
X = zeros(n,N);
X(:,1) = X0;
Y = zeros(ny,N);
U = zeros(m,N);
CGIters1 = zeros(N,1);
CGIters2 = zeros(N,1);

% LDIPM settings
mu_f = 1e-4;
mu_0 = 1e6;
maxIter = 150;
maxCGIter = 15000;
CGTol = 1e-6;

% Set options
opts.mu_f = mu_f;
opts.mu_0 = mu_0;
opts.maxIter = maxIter;
opts.printFlag = 0;
opts.maxCGIter = maxCGIter;
opts.CGTol = CGTol;

% Set parameters
gamma = 0.9;
opts.gamma = gamma;
opts.kappa = 2;

%% Run simulation
v0 = zeros(size(A_QP,1),1);
% v0 = -6*ones(size(A_QP,1),1);

for i = 1:N
     % Get current state vector
    Xi    = X(:,i);

    % Execute control law at current time-step - If unspecified flag is
    % passed in, the dynamics will integrate with no control law
    ref_i = 5;
    theta_i = [Xi; ref_i];
    c_QP = (W_MPC*theta_i);
    b_QP = b_MPC - L_MPC*theta_i;
   
    % * Get centered point for initialization
    % dTol =  1e-4;
    % [x_c,v_c,d] = logInteriorPoint_getCenteredPoint(H_QP,c_QP,A_QP,b_QP,v0,mu_0,1000,dTol);
    % v_init = v_c*0;
    v_init = v0;

    % * Run both LDIPM algorithms
    [x1,output1] = cgLDIPM_medstep(H_QP,c_QP,A_QP,b_QP,v_init,opts);
    [x2,output2] = cgLDIPM_medstep_inexact(H_QP,c_QP,A_QP,b_QP,v_init,opts);
    CGIters1(i) = sum(output1.CGIters);
    CGIters2(i) = sum(output2.CGIters);

    % Integrate dynamics with combined controls
    Ui =  x1(1:m);
    if i < N
        X(:,i+1) = A*Xi + B*Ui;
    end
    U(:,i)   = Ui;
    Y(:,i) = CPrime*Xi + DPrime*Ui;

    % Save the first data set for plotting later
    if i == 1  
        sampleData.H = H_QP;
        sampleData.c = c_QP;
        sampleData.A = A_QP;
        sampleData.b = b_QP;
        sampleData.U = x1;
        sampleData.v_init = v_init;
        sampleData.opts = opts;
        sampleData.output1 = output1;
        sampleData.output2 = output2;
    end
end

save('MPC_Data')

