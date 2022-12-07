%% prepare trajectory dataset
% N > (n_int+1)*T-1 = 19
%
% Yu Zheng, Florida State University
% 12/06/2022

clear all
clc

%% 1. load system matrices
T_sample = 0.01;   % sample time step
A_bar_d = [0.9950 0.0998;
           -0.0998 0.9950];
B_bar_d = [0.1;
           0.1];
[n_states,n_int] = size(B_bar_d);

n_meas = 6;
% C_temp = (1/sqrt(n_states))*randn(n_meas,n_states);
% C_obsv_d = zeros(n_meas,n_states);
% for C_col = 1:n_states
%     C_obsv_d(:,C_col) = C_temp(:,C_col)/norm(C_temp(:,C_col));
% end
% save C_obsv_d.mat C_obsv_d
load('C_obsv_d.mat');

D_obsv_d = zeros(n_meas,n_int);

% some horizon params
x0 = [7;7];
T = 5;
N = 3*(n_int+1)*T;


%% 2. runing parameters
N_samples      = 800;                % The total number of samples to run
T_final        = N_samples*T_sample;  % Total time for simulation (20s)


%% 3. run sim
u = -5 + 10*rand(N_samples,1);
u_traj = u(1:N);
HL_u = Get_Hanker(u_traj,T);

u_time = [linspace(0,T_final,N_samples).', u];
out = sim("sample_sys.slx");

y_traj_full  = out.logsout.getElement('y').Values.Data; 
y_traj = y_traj_full(1:N,:);
HL_y = Get_Hanker(y_traj,T);

%% 4.save trajectory data
save('traj.mat','u_traj','y_traj','-v7.3');
