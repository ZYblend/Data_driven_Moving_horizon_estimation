%% Runing model
%
% Yu Zheng, RASLab, FAMU-FSU College of Engineering, Tallahassee, 2021, Aug.

clear all
clc

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
C = [1 0];


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
U0          = zeros(n_int,L);


%% 4. runing parameters
tot_samples      = 800;               % The total number of samples to run
T_final        = T_sample*tot_samples;  % Total time for simulation (20s)
T_start_opt = T_sample*(L+1);           % start time for MHE
N_samples = 30;

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

%% Full information estimation

%% 6. data-driven MHE
rho = 0.95;
P = 10;
R = 1;

[H,f,Aeq,A_ineq,b_ineq] = Get_MHE_Param(P,R,rho,L,n_states,n_meas,u_d,x_d,y_d);

% n_alpha = size(HL_u,2);
% 
% % Error on y
% Aeq = [HL_u zeros(size(HL_u,1),n_states*L) zeros(size(HL_u,1),n_meas*L);
%        HL_y zeros(size(HL_y,1),n_states*L) eye(n_meas*L);
%        HL_x -eye(n_states*L) zeros(n_states*L,n_meas*L)];
% 
% H1 = zeros(1,n_alpha+n_states*L + n_meas*L);
% 
% H21 = [zeros(n_states,n_alpha) (sqrt(rho)^L)*sqrt(P)*eye(n_states) zeros(n_states,n_states*(L-1)) zeros(n_states,n_meas*L)];
% H22 = zeros(n_states*(L-1),n_alpha+n_states*L+ n_meas*L);
% H2 = [H21; H22];
% 
% rho_forget = zeros(n_meas*L);
% for idx = 1:L
%     rho_forget(n_meas*(idx-1)+1:n_meas*idx,n_meas*(idx-1)+1:n_meas*idx) = sqrt(rho)^(L+1-idx)*eye(n_meas);
% end
% H3 = [zeros(n_meas*L,n_alpha+n_states*L) rho_forget*sqrt(R)];
% 
% H = [H1;
%      H2;
%      H3];
% 
% f = [zeros(1,n_states);
%      eye(n_states);
%      zeros(n_states*(L-1),n_states);
%      zeros(n_meas*L,n_states)];
% 
% A_ineq = -[zeros(1,n_alpha) zeros(1,n_states*L) zeros(1,n_meas*L);
%           zeros(n_states*L,n_alpha) eye(n_states*L) zeros(n_states*L,n_meas*L);
%           zeros(n_meas*L,n_alpha) zeros(n_meas*L,n_states*L) zeros(n_meas*L)];
% b_ineq = zeros(1+n_states*L+n_meas*L,1);

% % Error on x
% Aeq = [HL_u zeros(size(HL_u,1),n_states*L) zeros(size(HL_u,1),n_states*L);
%        HL_y zeros(size(HL_y,1),n_states*L) zeros(n_meas*L,n_states*L);
%        HL_x -eye(n_states*L) eye(n_states*L)];
% 
% H1 = zeros(1,n_alpha+n_states*L + n_states*L);
% 
% H21 = [zeros(n_states,n_alpha) (rho^L)*sqrt(P)*eye(n_states) zeros(n_states,n_states*(L-1)) zeros(n_states,n_states*L)];
% H22 = zeros(n_states*(L-1),n_alpha+n_states*L+ n_states*L);
% H2 = [H21; H22];
% 
% rho_forget = zeros(n_states*L);
% for idx = 1:L
%     rho_forget(n_states*(idx-1)+1:n_states*idx,n_states*(idx-1)+1:n_states*idx) = rho^(L+1-idx)*eye(n_states);
% end
% H3 = [zeros(n_states*L,n_alpha+n_states*L) rho_forget*sqrt(R)];
% 
% H = [H1;
%      H2;
%      H3];
% 
% f = [zeros(1,n_states);
%      eye(n_states);
%      zeros(n_states*(L-1),n_states);
%      zeros(n_states*L,n_states)];
% 
% A_ineq = -[zeros(1,n_alpha) zeros(1,n_states*L) zeros(1,n_states*L);
%           zeros(n_states*L,n_alpha) eye(n_states*L) zeros(n_states*L,n_states*L);
%           zeros(n_states*L,n_alpha) zeros(n_states*L,n_states*L) zeros(n_states*L)];
% b_ineq = zeros(1+n_states*L+n_states*L,1);

