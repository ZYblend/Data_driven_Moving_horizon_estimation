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
A_bar_d = readmatrix('A_bar_d.csv');
B_bar_d = readmatrix('B_bar_d.csv');
C_obsv_d = readmatrix('C_obsv_d.csv');
D_obsv_d = readmatrix('D_obsv_d.csv');

[n_states,n_int] = size(B_bar_d);
n_meas = size(C_obsv_d,1);

% critical measurement matrix
Cm = ones(1,n_states);
n_critical = size(Cm,1);


%% 2. Check observability and controllabiltiy
disp('controllability')
disp(rank(ctrb(A_bar_d,B_bar_d))) % fully controllable with PID controller
disp('observability')
disp(rank(obsv(A_bar_d,C_obsv_d))) % fully observable


%% 3. Get T-horizon system matrices
T = 2*n_states;
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


%% 5. Pole placement for gain Control design
Pc = linspace(0.1,0.5, n_states);
K = place(A_bar_d,B_bar_d,Pc);
disp('discrete controller (A-B*K) eigenvalues: less than 1?')
disp(eig(A_bar_d-B_bar_d*K).')


%% 6. Simulation Initialization
% state initialization
x0          = 0.3+ 0.2*rand(n_states,1);
x0_hat      = zeros(n_states,1);
xd          = 0.5*ones(n_states,1);

yc_d = Cm*xd;

% Delay tape initialization
X0          = zeros(n_states,T);
Y0          = zeros(n_meas,T);
U0          = zeros(n_int,T+1);

offset = -inv(B_bar_d)*(A_bar_d-eye(n_states))*xd;


%% 7. runing parameters
N_samples      = 800;                % The total number of samples to run
T_final        = N_samples*T_sample;  % Total time for simulation (20s)
T_start_attack = 0.1*T_final;         % start injecting attack at 10s
T_start_opt(:)    = 1.5*T*T_sample;   % start state estimation at 
T_stop_attack = T_final;        % stop injecting attack at 10s

%% 8. Attack Parameters
T_start_attack = .2*T_final;  % Time to begin attack. Neede to sshow system responses with/without attacks in the same simulation
n_attack =  round(0.1*n_meas);
BDD_thresh = 5;  % Bad data detection tolerance

%% 9. moving-horizon FDIA parameters
I = sort(randperm(n_meas,n_attack)).';
% I = load('attack_support.mat').I;
I_attack_ini = repmat(I,1,T);           % fixed attack support
% flat attack support for T horizon
I_aux = linspace(0,(T-1)*n_meas,T);
I_aux = repmat(I_aux,n_attack,1);
I_attack = I_attack_ini+I_aux;
I_attack = I_attack(:);

% initial T attacks
e_ini = L2_FDIA(H0,I_attack(:),0.3*T*BDD_thresh);
e_ini_matrix = reshape(e_ini,n_meas,T);

% optimization solver parameters
[U,S,V] = svd(H0);
U_history = U(1:(T-1)*n_meas,:);
U_i = U((T-1)*n_meas+1:end,:);
S1 = S(1:n_states,1:n_states);

IO = [eye(n_states) zeros(n_states,T*n_meas-n_states)];
OI = [zeros(T*n_meas-n_states,n_states) eye(T*n_meas-n_states)];

A1 = inv(S1)*IO*U_i.';
A2 = OI*U_i.';
b1p = inv(S1)*IO*U_history.';
b2p = OI*U_history.';

max_iter = 5000;   % maximal number of iteration of PGA algorithm

%% 10. data-driven L2 MHE
% load trajectory data
traj = load('traj.mat');
u_traj = traj.u_traj;
y_traj = traj.y_traj;

HL_u = Get_Hanker(u_traj,T);
HL_y = Get_Hanker(y_traj,T);

% selection matrix
n_select = 2*n_states;
% E_idx = randperm(n_meas,n_select);
E_idx = load('y_selection.mat').E_idx;
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

n_alpha = size(HL_u,2);
Aeq = [zeros(size(HL_u,1),T*n_states), zeros(size(HL_u,1),n_meas*T), HL_u;
       zeros(size(HL_y,1),T*n_states), eye(n_meas*T), HL_y;
       CE_MH, E_MH, zeros(size(E_MH,1),n_alpha)];
Aeq_inv = pinv(Aeq,0.01);

%% 11. Data-driven L1 MHE
A_in1 = [zeros(size(HL_u,1),T*n_states), HL_u;
       zeros(size(HL_y,1),T*n_states), HL_y;
       CE_MH, zeros(size(E_MH,1),n_alpha)];
A_in2 = [zeros(size(HL_u,1),n_meas*T);
         eye(n_meas*T);
          E_MH];