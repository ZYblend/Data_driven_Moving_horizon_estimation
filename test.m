%% Test file for the relationship between 
%               Output reachability (markov matrix is full rank)
%               Persistent of excitation (Hanker matrix is full rank)
%
% Input is persistently exciting, and system is output reachable,
% Is the output trajectory is also persistently exciting?
%
% Yu Zheng, 2022/12/18
%

clear all
clc
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

D = zeros(2,2);

% A = [0.9950 0.0998;
%            -0.0998 0.9950];
% B = [0.1 -0.1;
%      -0.1   0.1];
% C = load('C_obsv_d.mat').C_obsv_d;
% 
% load('selection_index.mat');
% C = C(E_idx,:);
% C=rand(2,2);


n_states = size(A,1);
n_meas = size(C,1);
n_int = size(B,2);

D = zeros(n_meas,n_int);

% Check observability and controllabiltiy
disp('controllability')
disp(rank(ctrb(A,B))) % fully controllable 
disp('observability')
disp(rank(obsv(A,C))) % fully observable

T = 5;

% check if the system is output reachable
M = zeros(n_meas,(T+1)*n_int);
M (:,1:n_int) = D;
for idx_markov = 1:T
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
T_sample  = 0.01;
N_samples      = 30;                % The total number of samples to run
T_final        = N_samples*T_sample;  % Total time for simulation (20s)

u_traj = -10 + 20*rand(N_samples,n_int);
u_time = [linspace(0,T_final,N_samples).', u_traj];
x0 = rand(n_states,1);

out = sim("test_sys.slx");

y_traj  = out.logsout.getElement('y').Values.Data;


% test Hanker matrix
disp("input")
HL_u = Get_Hanker(u_traj,T+n_states);
disp("output")
HL_y = Get_Hanker(y_traj,T);