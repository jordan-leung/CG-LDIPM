clear
close all
clc

addpath('/Users/jordan/Documents/GitHub/optimizationFunctions');
addpath('/Users/jordan/Documents/GitHub/extraFunctions')
addpath('/Users/jordan/Documents/GitHub/terminalUnconstrainedGovernor');

set(groot,'defaultFigureCreateFcn',@(fig,~)addToolbarExplorationButtons(fig))

%% * Set simulation parameters

% Set run flag
% RUN_FLAG = 'standard';
% RUN_FLAG = 'standard_slack';
RUN_FLAG = 'fg';
% RUN_FLAG = 'fg_init';
% RUN_FLAG = 'ideal_fg';
% RUN_FLAG = 'ramp';

% Save figures and data?
saveFigFlag = 0;
saveDataFlag = 1;
dataFlag = 1;

% * Time vector
dt = 0.02;
t = 0:dt:10;
N = length(t);

% * Initial Conditions
X0 = [0; 0; 0; 0];
n = length(X0);

% * Define commanded reference
nr = 1;

% Case 1 reference
rConst =  1;
r =  rConst*ones(nr,N);


mCart  = 1;      % kg, mass of cart
mPend  = 0.1;      % kg, mass of pendulum
bf     = 0.1;      % N/m/s, damping/fraction coefficient
ell    = 1;      % m, length to pendulum centre of mass from bolt point
Im     = 1/3*mPend*ell^2;    % kg*m^2, mass moment of inertia of pendulum
g      = 9.81;     % m/s^2, gravity


% State Vector (X): x       = cart position (m) 
%                   xDot    = cart speed (m/s)
%                   phi     = pendulum angle pertubation (theta - pi) (rad)
%                   phiDot  = pendulum angular velocity (rad/s)


% Continuous-time dynamics
denom = Im*(mPend + mCart) + mPend*mCart*ell^2;

Ac = [0,                          1,                                 0,   0;...
     0, -(Im+mPend*ell^2)*bf/denom,             mPend^2*g*ell^2/denom,   0;...
     0,                          0,                                 0,   1;...
     0,        -mPend*ell*bf/denom,   mPend*g*ell*(mPend+mCart)/denom,   0];

 
Bc = [                     0;...
     (Im+mPend*ell^2)/denom;...
                          0;...
            mPend*ell/denom];
        
% Get discrete state-space matrices
m = size(Bc,2);
[A,B] = c2d(Ac,Bc,dt);



%% Generate tracking system matrices and constraint sets

% Min-max constraints
xmax = [Inf; Inf; 5*pi/180; Inf];
xmin = -xmax;
umax = 1;
umin = -umax;

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

% Set horizon
if strcmp(RUN_FLAG,'standard') || strcmp(RUN_FLAG,'standard_slack')
    N_MPC = 40;
    lambda = 1; % (141, 1000, 10 1 10 1), N < 44
    tRamp = 1; % placeholder
else
    N_MPC = 43; % 43
    lambda = 25; % 25, 50, 500
    tRamp = 2.3;  % 3.8 for (43,25), 2.9 for (43,50), 2.3 for (43,500)
end

T_MPC = dt;

% % * Set P, Q, R
Q = 1*eye(n);
Q(1,1) = 10;
Q(2,2) = 1;
Q(3,3) = 10;
Q(4,4) = 1;
R = 0.1; % Control effort cost

Q2 = Q;
R2 = R;

% % % Set Q and R for the terminal controller
% Q2 = Q;
% Q2(1,1) = 5;
% Q2(2,2) = 0.01;
% Q2(3,3) = 120;
% Q2(4,4) = 0.01;
% R2 = 0.5;
% Q = eye(n);
% R = 1;
% Q2 = Q;
% R2 = R;

K = dlqr(A,B,Q2,R2);
lambda_1 = 1;
P = dlyap((A-B*K)',lambda_1*(Q + K'*R*K));
P = (P+P')/2;

if mod(T_MPC,dt) ~= 0
    error('Mismatch in divisibility of T and dt')
end




%% Execute the simulation

% * Add controller arguments to a structure
controlArgs.n = n;
controlArgs.m = m;
controlArgs.A = A;
controlArgs.B = B;
controlArgs.P = P;
controlArgs.Q = Q;
controlArgs.R = R; 
controlArgs.K = K;
controlArgs.N = N_MPC; % Horzion lengt
controlArgs.T = T_MPC; % Sampling period
controlArgs.r = r;
controlArgs.xmin = xmin;
controlArgs.xmax = xmax;
controlArgs.umin = umin;
controlArgs.umax = umax;
controlArgs.Gx = Gx;
controlArgs.Gu = Gu;
controlArgs.Gr = Gr;
controlArgs.Yx = [0 0 1 0;...
      0 0 -1 0];
controlArgs.bx = [xmax(3); -xmin(3)];
controlArgs.Yu = [1; -1];
controlArgs.bu = [umax; - umin];
controlArgs.lambda = lambda;
controlArgs.MaxIter = 1e6;
controlArgs.XTol = 1e-8;
controlArgs.tRamp = tRamp;
% controlArgs.options = options;
