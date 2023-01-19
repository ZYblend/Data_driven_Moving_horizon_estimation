%% Runing model
%
% Yu Zheng, RASLab, FAMU-FSU College of Engineering, Tallahassee, 2023, Jan.

clear all
clc

%% 1. load system matrices
T_sample = 0.01;   % sample time step
A = [0.921 0 0.041 0;
     0 0.918 0 0.033;
     0 0 0.924 0;
     0 0 0 0.937];
B = [0.017 0.001;
           0.001 0.023;
           0 0.061;
           0.072 0];
C = [1 0 0 0;
     0 1 0 0];
% C = rand(6,4);
% A = [0.995 0.0998;
%      -0.0998 0.995];
% B = [0.1; 0.1];
% C = [1 0];


[n_states,n_int] = size(B);
n_meas = size(C,1);

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
X0_hat      = kron(ones(L,1),x0_hat);

% Delay tape initialization
X0          = zeros(n_states,L);
Y0          = zeros(n_meas,L);
U0          = zeros(n_int,L+1);


%% 4. runing parameters
tot_samples      = 800;                 % The total number of samples to run
T_final        = T_sample*tot_samples;  % Total time for simulation (20s)
T_start_opt = T_sample*(L+1);           % start time for MHE
N_samples = 26;                         % number of offline data

%% 5. Get trajectory data
u_traj = -5 + 10*rand(tot_samples,n_int);
u_time = [linspace(0,T_final,tot_samples).', u_traj];

out = sim("test_sys.slx");

y_traj  = out.logsout.getElement('y').Values.Data;
x_traj  = out.logsout.getElement('x').Values.Data;

% Historical data
u_d = u_traj(1:N_samples,:);
y_d = y_traj(1:N_samples,:);
x_d = x_traj(1:N_samples,:);

% Hankel matrix
HL_u = Get_Hankel(u_d,L);
HL_x = Get_Hankel(x_d,L);
HL_y = Get_Hankel(y_d,L);

C_bar = kron(eye(L),C);
Nc = null(C_bar);

HL_x_full = (eye(size(Nc,1))+ Nc*Nc.')*HL_x;

%% Check if HL(xd)*alpha = x
x_test1 = x_traj(N_samples+1:N_samples+L,:).';
y_test1 = y_traj(N_samples+1:N_samples+L,:).';
u_test1 = u_traj(N_samples+1:N_samples+L,:).';
alpha1 = [HL_u;HL_y;HL_x]\[u_test1(:);y_test1(:);x_test1(:)];

x_test2 = x_traj(N_samples+2:N_samples+L+1,:).';
y_test2 = y_traj(N_samples+2:N_samples+L+1,:).';
u_test2 = u_traj(N_samples+2:N_samples+L+1,:).';
alpha2 = [HL_u;HL_y;HL_x]\[u_test2(:);y_test2(:);x_test2(:)];


% HL_x*alpha - x_test(:)