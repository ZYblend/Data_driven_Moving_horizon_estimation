%% Runing model
%
% Yu Zheng, RASLab, FAMU-FSU College of Engineering, Tallahassee, 2021, Aug.

% clear all
% clc

%% 1. load system matrices
T_sample = 0.01;   % sample time step
% A = [0.921 0 0.041 0;
%      0 0.918 0 0.033;
%      0 0 0.924 0;
%      0 0 0 0.937];
% B = [0.017 0.001;
%            0.001 0.023;
%            0 0.061;
%            0.072 0];
% C = [1 0 0 0;
%      0 1 0 0];
A = [0.995 0.0998;
     -0.0998 0.995];
B = [0.1; 0.1];
% C = [1 0];




[n_states,n_int] = size(B);
% n_meas = size(C,1);

n_meas = 6;
% C = rand(n_meas,n_states);
C = load('C.mat').C;
D = zeros(n_meas,n_int);


%% 2. Check observability and controllabiltiy
disp('controllability')
disp(rank(ctrb(A,B))) % fully controllable with PID controller
disp('observability')
disp(rank(obsv(A,C))) % fully observable


%% 3. Simulation Initialization
% state initialization
L = 5;                                % horizon
x0          = 7*ones(n_states,1);
x0_hat      = [1; 2];
X0_hat      = zeros(L*n_states,1);
X0_hat(1:n_states) = x0_hat;

% Delay tape initialization
X0          = zeros(n_states,L);
Y0          = zeros(n_meas,L);
U0          = zeros(n_int,L);
U0_MD       = zeros(n_int,L+1);


%% 4. runing parameters
tot_samples      = 800;               % The total number of samples to run
T_final        = T_sample*tot_samples;  % Total time for simulation (20s)
T_start_opt = T_sample*(L+1);           % start time for MHE

N_samples = 30;     % number of offline data

%% 5. Get trajectory data
u_traj = -5 + 10*rand(tot_samples,n_int);
u_time = [linspace(0,T_final,tot_samples).', u_traj];

out = sim("test_sys.slx");

y_traj  = out.logsout.getElement('y').Values.Data;
x_traj  = out.logsout.getElement('x').Values.Data;
u_traj  = out.logsout.getElement('u').Values.Data;

% Historical data
u_d = u_traj(1:N_samples,:);
y_d = y_traj(1:N_samples,:);
x_d = x_traj(1:N_samples,:);
zd = {u_d, x_d, y_d};

% Hankel matrix
HL_u = Get_Hankel(u_d,L);
HL_x = Get_Hankel(x_d,L);
HL_y = Get_Hankel(y_d,L);

%% Check if HL(xd)*alpha = x
x_test = x_traj(N_samples+1:N_samples+L,:).';
y_test = y_traj(N_samples+1:N_samples+L,:).';
u_test = u_traj(N_samples+1:N_samples+L,:).';
alpha = [HL_u;HL_y]\[u_test(:);y_test(:)];
HL_x*alpha - x_test(:)

% update u_time for simulation
u_sim = 10*rand(tot_samples,n_int);
u_time_sim = [linspace(0,T_final,tot_samples).', u_sim];


%% 6. data-driven MHE
n_alpha = size(HL_u,2);
z_hat0 = zeros(n_alpha+n_states*L + n_meas*L,1);
% L2
rho = 1;
P = 1;
R = 10;

[H_DDL2,f_DDL2,Aeq_DDL2,A_ineq_DDL2,b_ineq_DDL2] = Get_DDL2_MHE_Param(P,R,rho,L,n_states,n_meas,u_d,x_d,y_d);

% L1
[H_L1_DDL1,H_L2_DDL1,f_DDL1,Aeq_DDL1,A_ineq_DDL1,b_ineq_DDL1] = Get_DDL1_MHE_Param(P,R,rho,L,n_states,n_meas,u_d,x_d,y_d);

%% 7. Attack Parameters
T_start_attack = .1*T_final;  % Time to begin attack. Neede to sshow system responses with/without attacks in the same simulation
T_stop_attack = T_final;             % stop injecting attack at 8s
% n_attack =  round(0.2*n_meas);
% I = randperm(n_meas,n_attack);
indicator = zeros(n_meas,1);
indicator(I) = 1;

%% 8. Model based estimators
% L2
[H_MDL2,f_MDL2,Aeq_MDL2,H1_MDL2,F_MDL2] = Get_MDL2_MHE_Param(P,R,rho,L,n_states,n_meas,A,B,C);

% L1 
[H_L1_MDL1,H_L2_MDL1,f_MDL1,Aeq_MDL1,H1_MDL1,F_MDL1] = Get_MDL1_MHE_Param(P,R,rho,L,n_states,n_meas,A,B,C);
