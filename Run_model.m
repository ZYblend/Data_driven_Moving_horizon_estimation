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


%% 2. Check observability and controllabiltiy
disp('controllability')
disp(rank(ctrb(A_bar_d,B_bar_d))) % fully controllable with PID controller
disp('observability')
disp(rank(obsv(A_bar_d,C_obsv_d))) % fully observable


%% 3. Get T-horizon system matrices
T = 5;
[H0,H1,F] = opti_params(A_bar_d,B_bar_d,C_obsv_d,T);
% H0: state-ouput linear map                       [n_meas*T-by-n_states]
% H1:  input-output linear map                     [n_meas*T-by-n_int*(T-1)]
% F:  Observer input-state propagation matrix      [n_meas-by-n_int*(T-1)]


%% 4. Get observer gains for L2 observer
H0_pinv = pinv(H0,0.001);
Ly = mpower(A_bar_d,T)*H0_pinv;
Lu = F-mpower(A_bar_d,T)*H0_pinv*H1;

% residual
H0_perp = eye(size(H0,1)) - H0*H0_pinv;


%% 5. Simulation Initialization
% state initialization
x0          = [7;7];
x0_hat      = [1;2];

% Delay tape initialization
X0          = zeros(n_states,T);
Y0          = zeros(n_meas,T);
U0          = zeros(n_int,T+1);


%% 6. runing parameters
N_samples      = 800;                % The total number of samples to run
T_final        = N_samples*T_sample;  % Total time for simulation (20s)
T_start_attack = 0.1*T_final;         % start injecting attack at 10s
T_start_opt(:)    = 1.5*T*T_sample;   % start state estimation at 
T_stop_attack = T_final;        % stop injecting attack at 10s

%% 7. inputs trajectory
% u = -10 + 20*rand(N_samples,1);
% u_time = [linspace(0,T_final,N_samples).', u];
load u_time.mat

% %% 8. Attack Parameters
T_start_attack = .2*T_final;  % Time to begin attack. Neede to sshow system responses with/without attacks in the same simulation
n_attack =  round(0.2*n_meas);
BDD_thresh = 5;  % Bad data detection tolerance
I = randperm(n_meas,n_attack);
indicator = zeros(n_meas,1);
indicator(I) = 1;

% 
%% 10. data-driven MHE
% load trajectory data
traj = load('traj.mat');
u_traj = traj.u_traj;
y_traj = traj.y_traj;

HL_u = Get_Hanker(u_traj,T);
HL_y = Get_Hanker(y_traj,T);

% selection matrix
% n_select = 2*n_states;
% E_idx = randperm(n_meas,n_select);
load('selection_index.mat');
% E_idx = [1,2,3,4,5,6];
n_select = length(E_idx);
E = zeros(n_select,n_meas);
for iter = 1:n_select
    E(iter,E_idx(iter)) = 1;
end
E_MH = kron(eye(T),E);

% model for the selected measurements
C_E = C_obsv_d(E_idx,:);
D_E = D_obsv_d(E_idx,:);
disp('observability of selected measurements')
disp(rank(obsv(A_bar_d,C_E))) % fully observable
CE_MH = kron(eye(T),C_E);
DE_MH = kron(eye(T),D_E);

%% L2
n_alpha = size(HL_u,2);
H = [zeros(size(HL_u,1),T*n_states), HL_u;
       zeros(size(HL_y,1),T*n_states), HL_y;
       CE_MH, zeros(size(E_MH,1),n_alpha)];
H_pinv = pinv(H,0.001);

%% L1 MHE
A_in1 = H;
A_in2 = [zeros(size(HL_u,1),n_meas*T);
         eye(n_meas*T);
         E_MH];