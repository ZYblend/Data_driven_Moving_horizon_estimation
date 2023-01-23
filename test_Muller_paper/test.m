%% Test file for the relationship between 
%               Output reachability (markov matrix is full rank)
%               Persistent of excitation (Hanker matrix is full rank)
%
% Input is persistently exciting, and system is output reachable,
% Is the output trajectory is also persistently exciting?
%
% Yu Zheng, 2022/12/18
%

clear 
clc

%% System Dynamics
% A = [0.921 0 0.041 0;
%      0 0.918 0 0.033;
%      0 0 0.924 0;
%      0 0 0 0.937];
% B = [0.017 0.001;
%      0.001 0.023;
%      0 0.061;
%      0.072 0];
% C = [1 0 0 0;
%      0 1 0 0];

% D = zeros(2,2);

% A = [0.9950 0.0998;
%            -0.0998 0.9950];
% B = [0.1 -0.1;
%      -0.1   0.1];
% C = load('C_obsv_d.mat').C_obsv_d;
% 
% load('selection_index.mat');
% C = C(E_idx,:);
% C=rand(2,2);

A = [0.4 -0.3 0 0.1;
     -0.3 0 0.8 -0.1;
     0.1 -0.7 -0.4 0;
     0.2 -0.5 0.5 0.4];

B = [0; -1; 1.4; 0];
C = [-0.7 0 -2 0.4;
     0.5 0.6 0.7 0];
D = 0.2;


 
n_states = size(A,1);
n_meas = size(C,1);
n_int = size(B,2);

% D = zeros(n_meas,n_int);

% Check observability and controllabiltiy
disp('controllability')
disp(rank(ctrb(A,B))) % fully controllable 
disp('observability')
disp(rank(obsv(A,C))) % fully observable

L = 2; % Horizon lenth

% check if the system is output reachable
M = zeros(n_meas,(L+1)*n_int);
M (:,1:n_int) = D;
for idx_markov = 1:n_states
    M(:,n_int*(idx_markov)+1:n_int*(idx_markov+1)) = C*mpower(A,idx_markov-1)*B;
end
disp('rank of Markov matrix')
disp(rank(M))
if rank(M)==size(M,1)
    disp('the system is output reachable')
else
    disp('the system is not output reachable')
end

% get trajetory with persistently exciting inputs
T_sample       = 0.01;
N_samples      = 24;                % The total number of samples to run
N_add_samples  = 10;                % Number of additional samples after initial runs
N_total_sample = N_samples + N_add_samples;
T_final        = (N_total_sample)*T_sample;  % Total time for simulation (20s)

u_traj_inp = -10 + 20*rand(N_total_sample,n_int);
u_time = [linspace(0,T_final,N_total_sample).', u_traj_inp];
x0 = rand(n_states,1);

out = sim("test_sys.slx");

y_traj  = out.logsout.getElement('y').Values.Data;
x_traj  = out.logsout.getElement('x').Values.Data;
u_traj  = out.logsout.getElement('u').Values.Data;


% Historical data
u_d = u_traj(1:N_samples,:);
y_d = y_traj(1:N_samples,:);
x_d = x_traj(1:N_samples,:);

%% Testing Hanker matrix adjoint property
disp("*************************************")
disp('HANKEL MATRIX ADJOINT PROPERTY')
disp("*************************************")

disp("input:")
HL_u = Get_Hankel(u_d,L);
alpha_u = rand(size(HL_u,2),1);
u = HL_u*alpha_u;
u_d_T = u_d.';
HL_u_adj = Get_Hankel_adj(alpha_u,L,n_int);
disp([HL_u_adj*u_d_T(:) HL_u*alpha_u])
disp("Residual =" + num2str(norm(HL_u_adj*u_d_T(:)-HL_u*alpha_u)))
disp("*************************************")

disp("*************************************")
disp("output")
HL_y = Get_Hankel(y_d,L);
alpha_y = rand(size(HL_y,2),1);
y = HL_y*alpha_y;
y_d_T = y_d.';
HL_y_adj = Get_Hankel_adj(alpha_y,L,n_meas);
disp([HL_y_adj*y_d_T(:) HL_y*alpha_y])
disp("Residual =" + num2str(norm(HL_y_adj*y_d_T(:)-HL_y*alpha_y)))
disp("*************************************")

disp("*************************************")
disp("state")
HL_x = Get_Hankel(x_d,L);
alpha_x = rand(size(HL_x,2),1);
x = HL_x*alpha_x;
x_d_T = x_d.';
HL_x_adj = Get_Hankel_adj(alpha_x,L,n_states);
disp([HL_x_adj*x_d_T(:) HL_x*alpha_x])
disp("Residual =" + num2str(norm(HL_x_adj*x_d_T(:)-HL_x*alpha_x)))
disp("*************************************")


t_start = N_samples+1;  % start time for windowed mesaurements
alphas = zeros(size(HL_u,2),N_total_sample-t_start+1);
alphas2 = zeros(size(HL_u,2),N_total_sample-t_start+1);
% x_pred = zeros(n_states,N_total_sample-t_start+1);  % predicted states from Hankel matrices

X_win = zeros(L*n_states,N_total_sample-t_start+1);
Y_win = zeros(L*n_meas,N_total_sample-t_start+1);
for iter = t_start:N_total_sample
    curr_win = iter-L+1:iter;  % current window

    % input, output and states at current window
    u_curr = u_traj(curr_win,:).';
    y_curr = y_traj(curr_win,:).';
    x_curr = x_traj(curr_win,:).';

    X_win(:,iter-t_start+1) = x_curr(:);
    Y_win(:,iter-t_start+1) = y_curr(:);
    
    % Model verification
    alphas(:,iter-t_start+1) = [HL_u;HL_y]\[u_curr(:);y_curr(:)];
    alphas2(:,iter-t_start+1) = [HL_u;HL_y;HL_x]\[u_curr(:);y_curr(:);x_curr(:)];

end
    
% Verify Hankel
C_bar = kron(eye(L),C);
D_bar = kron(eye(L),D);
% Y_win - C_bar*X_win
% Y_win - HL_y*alphas
X_win - HL_x*alphas
C_bar*(X_win - HL_x*alphas)


