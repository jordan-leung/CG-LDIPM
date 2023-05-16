clear
close all
clc

addpath('/Users/jordan/Documents/GitHub/MATLAB/optimizationFunctions');
addpath('/Users/jordan/Documents/GitHub/MATLAB/CG-LDIPM/InexactNewton')

set(groot,'defaultFigureCreateFcn',@(fig,~)addToolbarExplorationButtons(fig))

%% * Set simulation parameters

% Set horizon
N_MPC = 15; % 104 is the minimum for X0 = 0, r = 3
T_MPC = 0.1;

MASS = 1724; % vehicle mass, kg
Izz = 1100; % yaw moment-of-intertia, kg/m^2
a_cg = 1.35; % front axle-CG distance, m
b_cg = 1.15; % rear axle-CG distance, m
C_alphaf = 90000; % front corning stiffness, N/rad
C_alphar = 138000; % front corning stiffness, N/rad
Ux = 10; % longitudinal velocity, m/s


% State Vector (X): beta       = velocity ratio, Uy/Ux
%                   r          = yaw rate (rad/s)
%                   y          = y position (m)

% Control Vector:   u          = steering angle (rad)

% * Initial Conditions
X0 = [0; 0; 0];
n = length(X0);

% * Define commanded reference
nr = 1;
r =  0;

% Continuous-time dynamics

A = [ -(C_alphaf + C_alphar)/(MASS*Ux),          -(a_cg*C_alphaf-b_cg*C_alphar)/(MASS*Ux^2)-1,  0;...
      -(a_cg*C_alphaf - b_cg*C_alphar)/(Izz),    -(a_cg^2*C_alphaf + b_cg^2*C_alphar)/(Izz*Ux), 0;
      Ux,                                        0                                              0];

B = [ C_alphaf/(MASS*Ux);...
      a_cg*C_alphaf/Izz;...
      0];
          
% Get discrete state-space matrices
m = size(B,2);
[Ad,Bd] = c2d(A,B,T_MPC);



%% * Set control parameters and run the simulation


% * Set P, Q, R
Q = diag([1,1,10]); %  State error cost
R = 1; % Control effort cost
[K,P] = dlqr(Ad,Bd,Q,R);

% * Set constraints 
umax = 1;
umin = -1;
xmax = [0.2; 4; 4]; 
% xmax = [0.05; 4; 4]; % 0.7854, 3.4907
xmin = -xmax;
MaxIter = 10000; % 7 and 6200 at 0.2 s, R = 0.1, 10 and 8000 at 0.15 s, 15 and 5500 at 0.1 s
xTol_MPC = 1e-8;  % ||x||_inf optimization tolerance for MPC optimizer


%% Generate Onfinity set....

% % constraints in the form Cx + Du \in Y

CPrime = [eye(n);zeros(m,n); -eye(n); zeros(m,n)];
DPrime = [zeros(n,m); eye(m); zeros(n,m); -eye(m)];

% Tracking sttates in the form z = Ex + Fu
Es = [0 0 1];
F = 0;

[nx,nu] = size(Bd);

Z = [eye(nx)-Ad, Bd, zeros(nx,nu); Es, F, -eye(nr)]; 
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

%++++++++++++++++++++++++ Generate Oinf +++++++++++++++++++++++++++++++++++
A_bar = Ad-Bd*K;
B_bar = Bd*K*Gx + Bd*Gu;
C_bar = CPrime-DPrime*K;
D_bar = DPrime*K*Gx + DPrime*Gu;
  
Ay = eye(2*(n+m));
by = [xmax; umax; -xmin; -umin];

% Generate Oinf and flips column order from [r x] to [x r]
[T_inf, c_inf] = genOinf(A_bar, B_bar, C_bar, D_bar , Ay, by, 200,.01,1e-12); % generate Oinf
temp = T_inf(:,nr+1:end);
T_inf = [temp, T_inf(:,1:nr)];
% save('Oinf','T_inf','c_inf')
% load('Oinf')

% Plot the slice of Gamma_N(r) in the (beta,y) space at r = 0

% Generate the constraint set Mz + Lx \leq b
POinf=Polyhedron('A',T_inf,'b',c_inf);
Tx = T_inf(:,1:n);
Tr = T_inf(:,n+1:end);
r0 = r(1);
b_inf = c_inf - Tr*r0;
Oinf = Polyhedron('A',Tx,'b',b_inf);
Y = Polyhedron('A',Ay,'b',by);
[M,L,b] = feas_set_matrices(N_MPC,Ad,Bd,CPrime,DPrime,Y,Oinf);
nn = size(L,2);
mm = size(M,2);

% * FOR PLOTTING A SINGLE REGION
% perform the projection
Pi = [eye(nn),zeros(nn,mm)];
H = polyh(struct('B',[L M],'b',b),'h');
PP = im(H,Pi);
PP = hrep(PP);
Gamma_A = PP.B;
Gamma_b = PP.b;
Gamma = Polyhedron('A',PP.B,'b',PP.b);
box on; grid on;

% Get the slice
rFix = 0;
GammaSlice_A = Gamma_A(:,[1 3]);
GammaSlice_b = Gamma_b - Gamma_A(:,2)*rFix; %  redundant i suppose
GammaSlice = Polyhedron('A',GammaSlice_A,'b',GammaSlice_b);

% Also get the slice of O_infty
OinfSlice_A = Tx(:,[1 3]);
OinfSlice_b = b_inf - Tx(:,2)*rFix;
OinfSlice = Polyhedron('A',OinfSlice_A,'b',OinfSlice_b);


%% Simulation 


% * Add controller arguments to a structure
controlArgs.n = n;
controlArgs.m = m;
controlArgs.A = Ad;
controlArgs.B = Bd;
controlArgs.P = P;
controlArgs.Q = Q;
controlArgs.R = R; 
controlArgs.K = K;
controlArgs.N = N_MPC; % Horzion lengt
controlArgs.T = T_MPC; % Sampling period
controlArgs.MaxIter = MaxIter;
controlArgs.xTol = xTol_MPC;
controlArgs.umax = umax; % Control constraints
controlArgs.umin = umin;
controlArgs.T_inf = T_inf;
controlArgs.c_inf = c_inf;
controlArgs.r = r;
controlArgs.CPrime = CPrime;
controlArgs.DPrime = DPrime;
controlArgs.Ay = Ay;
controlArgs.by = by;
controlArgs.Gx = Gx;
controlArgs.Gu = Gu;
controlArgs.Gr = Gr;

% Initial reference command
controlArgs.nu_0 = X0(3); 

% LDIPM options
LDIPM_opts.mu_f = 1e-8;
LDIPM_opts.mu_0 = 1e8;
LDIPM_opts.mu_weight = 1; 
LDIPM_opts.mu_bound = [1e-10, 1e-2];
LDIPM_opts.maxIter = 1000;
LDIPM_opts.printFlag = 0;
LDIPM_opts.FNormLimit = 10; % bound on the neighborhood
controlArgs.LDIPM_opts = LDIPM_opts;

%% Running simple experiments

% Generate the constant QP matrices (fixed horizon lengths)
[H_MPC,W_MPC,M_MPC,L_MPC,b_MPC,AHat,BHat] = generateQPMatrices_tracking(N_MPC,Ad,Bd,CPrime,DPrime,P,Q,R,Ay,by,T_inf,c_inf,Gx,Gu);
H_QP = H_MPC;
A_QP = M_MPC;

% Loop over initial conditions. Fix r = 0 and just look over side-slip
% angle and angular velocity
NSample = 50;
betaSample = linspace(0.95*xmin(1),0.95*xmax(1),NSample);
ySample = linspace(0.95*xmin(3),0.95*xmax(3),NSample);
[BETA,YY] = meshgrid(betaSample,ySample);
NUM_ITER = BETA*0; % initialize
v0 = zeros(size(M_MPC,1),1);
for i = 1:NSample
    i
    for j = 1:NSample   
        % Define initial state and parameter
        X_i = [betaSample(i); 0 ; ySample(j)];
        theta_i = [X_i; r];
        
        % Check feasibility
        feasVec = Gamma_A*X_i - Gamma_b;
        feasFlag = 1;
        for k = 1:length(feasVec)
            if feasVec(k) > 1e-6
                feasFlag = 0;
                break
            end
        end
        
        % Run MPC if feasibile
        if feasFlag
            f_QP = (W_MPC*theta_i);
            b_QP = b_MPC - L_MPC*theta_i;
            [~,output] = logInteriorPoint_conjgrad_INB_edit(H_QP,f_QP,A_QP,b_QP,v0,LDIPM_opts);
            NUM_ITER(i,j) = sum(output.CGIters);
        else
            NUM_ITER(i,j) = NaN;
        end
    end
end

save('data_10p0')



