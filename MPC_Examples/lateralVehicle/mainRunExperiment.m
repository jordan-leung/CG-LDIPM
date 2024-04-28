clear
close all
clc

addpath('/Users/jordan/Documents/GitHub/MATLAB/optimizationFunctions');
addpath('/Users/jordan/Documents/GitHub/MATLAB/extraFunctions')
addpath('/Users/jordan/Documents/GitHub/CG-LDIPM/MPC_Examples')

set(groot,'defaultFigureCreateFcn',@(fig,~)addToolbarExplorationButtons(fig))

%% * Set simulation parameters


% Print flag
printFlag = 0;

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
X0 = [4; 0; 0; 0];
% X0 = [-4.9; 0; 0; 0];
n = length(X0);

% * Define commanded reference
nr = 1;

% Case 1 reference
r =  5*ones(nr,N);


% Continuous-time dynamics
A = [0, Vx, Vx, 0;...
     0, 0, 0, 1;...
     0, 0, -2*C_alpha/(MASS*Vx), C_alpha*(ell_r - ell_f)/(MASS*Vx^2) - 1;...
     0, 0, C_alpha*(ell_r - ell_f)/Izz, -C_alpha*(ell_r^2 + ell_f^2)/(Izz*Vx)];
 
 B = [0; 0; C_alpha/(MASS*Vx); C_alpha*ell_f/Izz];

% Get discrete state-space matrices
m = size(B,2);
[Ad,Bd] = c2d(A,B,dt);
T_MPC = dt;


%% Generate tracking system matrices and constraint sets

% % constraints in the form Cx + Du \in Y... Relates states to slip angles
CPrime = [0, 0, -1, -ell_f/Vx;...
          0, 0, -1, ell_r/Vx;...
          0, 0, 0, 0];
DPrime = [1; 0; 1];
yPrime = [8; 8; 30]*pi/180; 
EPrime = CPrime;
gPrime = yPrime;

%% Simulation and MPC parameters

% Set horizon
N_MPC = 30;

% * Set P, Q, R
Q = 0.1*eye(n);
Q(1,1) = 1;
R = 0.1; % Control effort cost
P = Q;


if mod(T_MPC,dt) ~= 0
    error('Mismatch in divisibility of T and dt')
end

  % LDIPM options
LDIPM_opts.mu_f = 1e-8;
LDIPM_opts.mu_0 = 1e8;
LDIPM_opts.mu_weight = 1; 
LDIPM_opts.dLim = 1 - 0.01; 
LDIPM_opts.mu_bound = [1e-10, 1e-2];
LDIPM_opts.kappa_scale = 1;
LDIPM_opts.maxIter = 150;
LDIPM_opts.printFlag = printFlag;
controlArgs.LDIPM_opts = LDIPM_opts;

% Generate Matrices
[H,F,M,L,b,AHat,BHat] = generateMatrices_forCGLDIPM(N_MPC,Ad,Bd,CPrime,DPrime,yPrime,P,Q,R,EPrime,gPrime);

%% Solve











