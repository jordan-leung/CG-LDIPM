clear
close all
clc

addpath('/Users/jordan/Documents/GitHub/optimizationFunctions');
addpath('/Users/jordan/Documents/GitHub/extraFunctions')

set(groot,'defaultFigureCreateFcn',@(fig,~)addToolbarExplorationButtons(fig))

%% * Set simulation parameters

% * Time vector
dt = 0.01;


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
theta_i = [X0; 5];
H_QP = H_MPC;
A_QP = M_MPC;
c_QP = (W_MPC*theta_i);
b_QP = b_MPC - L_MPC*theta_i;

       
%% Test and plot trajectory
U = quadprog(H_QP,c_QP,A_QP,b_QP);

% Get X trajectories
X = zeros(n,N_MPC);
X(:,1) = X0;
Y = zeros(ny,N_MPC);
Y(:,1) = CPrime*X0;
for i = 1:N_MPC
    if i < N_MPC
        X(:,i+1) = A*X(:,i) + B*U(i);
    end
    Y(:,i) = CPrime*X(:,i) + DPrime*U(i);
end

% Plot: s, psi, beta, alpha_f, alpha_r, delta
close all

set(0,'defaultLineLineWidth',2)
set(0,'defaultAxesFontSize',12)
labelsize = 15;

saveDirectory = 'C:\Users\Jordan Leung\Documents\MATLAB\PhD\19-08-26_HCWToyProblem';


colorMatrix = [0 0 0; 228 26 28; 55 126 184; 77 175 74]./255;
greyColor = [111 111 111]/255;

%% Producing plots
figSize = [0 0 0.3 0.3];
subPlotGap = [0.13 0.10];
subPlotH = [0.1 0.05];
subPlotW = [0.1 0.05];

t = 0:dt:dt*(N_MPC-1);

figure
% s
set(gcf,'units','normalized','position',figSize)
subtightplot(2,3,1,subPlotGap,subPlotH,subPlotW)
hold on; box on; grid on
plot(t,X(1,:),'Color',colorMatrix(1,:))
plot(t,t*0+5,'color',colorMatrix(2,:),'linestyle','-.')
xlabel('Time','interpreter','Latex','Fontsize',labelsize)
ylabel('${s}$ [m]','interpreter','Latex','FontSize',labelsize)
% ylim([-Inf 5.1])

% psi
subtightplot(2,3,2,subPlotGap,subPlotH,subPlotW)
hold on; box on; grid on;
plot(t,X(2,:)*180/pi,'Color',colorMatrix(1,:))
xlabel('Time','interpreter','Latex','Fontsize',labelsize)
ylabel('${\psi}$ [deg]','interpreter','Latex','FontSize',labelsize)

% beta
subtightplot(2,3,3,subPlotGap,subPlotH,subPlotW)
hold on; box on; grid on;
plot(t,X(3,:)*180/pi,'Color',colorMatrix(1,:))
xlabel('Time','interpreter','Latex','Fontsize',labelsize)
ylabel('${\beta}$ [deg]','interpreter','Latex','FontSize',labelsize)

% alpha_f
subtightplot(2,3,4,subPlotGap,subPlotH,subPlotW)
hold on; box on; grid on;
plot(t,t*0+ymax(1)*180/pi,'Linestyle','-.','LineWidth',2,'Color',greyColor)
plot(t,t*0+ymin(1)*180/pi,'Linestyle','-.','LineWidth',2,'Color',greyColor)
plot(t,Y(1,:)*180/pi,'k','LineWidth',2)
xlim([0 max(t)])
xlabel('Time','interpreter','Latex','Fontsize',labelsize)
ylabel('$\alpha_f$ [deg]','interpreter','Latex','FontSize',labelsize)

% alpha_r
subtightplot(2,3,5,subPlotGap,subPlotH,subPlotW)
hold on; box on; grid on;
plot(t,t*0+ymax(2)*180/pi,'Linestyle','-.','LineWidth',2,'Color',greyColor)
plot(t,t*0+ymin(2)*180/pi,'Linestyle','-.','LineWidth',2,'Color',greyColor)
plot(t,Y(2,:)*180/pi,'k','LineWidth',2)
xlim([0 max(t)])
xlabel('Time','interpreter','Latex','Fontsize',labelsize)
ylabel('$\alpha_r$ [deg]','interpreter','Latex','FontSize',labelsize)

% delta
subtightplot(2,3,6,subPlotGap,subPlotH,subPlotW)
hold on; box on; grid on;
plot(t,t*0+ymax(3)*180/pi,'Linestyle','-.','LineWidth',2,'Color',greyColor)
plot(t,t*0+ymin(3)*180/pi,'Linestyle','-.','LineWidth',2,'Color',greyColor)
plot(t,U*180/pi,'Color',colorMatrix(1,:))
xlabel('Time','interpreter','Latex','Fontsize',labelsize)
ylabel('${\delta}$ [deg]','interpreter','Latex','FontSize',labelsize)
% ylim([1.1*ymin(3) 1.1*ymax(3)])


saveFigFlag = 0;
if saveFigFlag == 1
    filename = strcat('./Plots/singleManeuver_cg');
    saveas(gcf,filename,'epsc')
    save('OUTPUT_CG')
end




%% Rename and save
H = H_QP;
A = A_QP;
c = c_QP;
b = b_QP;
save('MPC_Data','H','c','A','b');










