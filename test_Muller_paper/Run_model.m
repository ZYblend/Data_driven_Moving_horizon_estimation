%% Runing model
% Content:
%         1. load system matrices
%         2. Check observability and controllabiltiy
%         3. Get T-horizon system matrices
%         4. Get observer gains for L2 observer
%         5. Pole placement for gain Control design
%         6. Simulation Initialization
%         7. Runing parameters (time frame...)
%
% Yu Zheng, RASLab, FAMU-FSU College of Engineering, Tallahassee, 2021, Aug.

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
L = 7;                                % horizon
x0          = 7*eye(n_states,1);
x0_hat      = [1 2 1 2];

% Delay tape initialization
X0          = zeros(n_states,L);
Y0          = zeros(n_meas,L);
U0          = zeros(n_int,L+1);


%% 4. runing parameters
tot_samples      = 800;               % The total number of samples to run
T_final        = T_sample*tot_samples;  % Total time for simulation (20s)
T_start_opt = T_sample*(L+1);           % start time for MHE
N_samples = 100;

%% 5. Get trajectory data
u_traj = 20*rand(N_samples,n_int);
u_time = [linspace(0,T_final,N_samples).', u_traj];

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



%% 6. data-driven MHE
n_alpha = size(HL_u,2);
Aeq = [HL_u zeros(size(HL_u,1),n_states*L);
     HL_y zeros(size(HL_y,1),n_states*L);
     HL_x -eye(n_states*L)];
H1 = zeros(1,n_alpha+n_states*L);
H2 = [zeros(n_states,n_alpha) eye(n_states) zeros(n_states,n_states*(L-1))];
H3 = zeros(n_states*(L-1),n_alpha+n_states*L);
H = [H1;
     H2;
     H3];

f = [zeros(1,n_states);
     eye(n_states);
     zeros(n_states*(L-1),n_states)];
A_ineq = [zeros(1,n_alpha) zeros(1,n_states*L);
     zeros(n_states*L,n_alpha) eye(n_states*L)];
b_ineq = zeros(1+n_states*L,1);
