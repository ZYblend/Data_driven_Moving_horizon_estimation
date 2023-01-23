%% Numerican experiement for "A trajectory-based framework for data-driven system analysis and control"
%
% Yu Zheng, Florida State University
% yz19b@fsu.edu
% 01/19/2023

%% 1. system dynamics
A = [0.4 -0.3 0 0.1;
     -0.3 0 0.8 -0.1;
     0.1 -0.7 -0.4 0;
     0.2 -0.5 0.5 0.4];

B = [0; -1; 1.4; 0];
C = [-0.7 0 -2 0.4];
D = 0.2;

% check controllability and obervability
disp('controllability')
disp(rank(ctrb(A,B))) % fully controllable
disp('observability')
disp(rank(obsv(A,C))) % fully observable

% dimentions
[n_states,n_int] = size(B);
n_meas = size(C,1);


%% 2. System simulation parameters
T_sample = 0.01;   % sample time step
L  = 5;            % horizon
x0 = zeros(n_states,1);  % initial condition

tot_samples = 2000;               % The total number of samples to run
T_final     = T_sample*tot_samples;  % Total time for simulation (20s)
N_samples   = 30;                  % number of sample collected offline

%% 3. Get trajectory data
u_traj = -30 + 60*rand(tot_samples,n_int);
u_time = [linspace(0,T_final,tot_samples).', u_traj];

out = sim("test_sys.slx");

y_traj  = out.logsout.getElement('y').Values.Data;
x_traj  = out.logsout.getElement('x').Values.Data;

% Historical data
u_d = u_traj(1:N_samples,:);
y_d = y_traj(1:N_samples,:);
x_d = x_traj(1:N_samples,:);

% Hankel matrix
disp("input trajectory is");
HL_u = Get_Hankel(u_d,L);

disp(' ');
disp("state trajectory is");
HL_x = Get_Hankel(x_d,L);

disp(' ');
disp("output trajectory is");
HL_y = Get_Hankel(y_d,L);

%% 4. Check if HL(xd)*alpha = x
x_test = x_traj(N_samples+1:N_samples+L,:).';
y_test = y_traj(N_samples+1:N_samples+L,:).';
u_test = u_traj(N_samples+1:N_samples+L,:).';

% by william fundamental lemma
alpha = [HL_u;HL_y]\[u_test(:);y_test(:)];

disp(' ');
disp("HL(x)*\alpha - x = ")
disp(num2str( HL_x*alpha - x_test(:) ));

disp(' ');
CL = kron(eye(L),C);
disp("C(HL(x)*\alpha - x) =")
disp(num2str( CL*(HL_x*alpha - x_test(:)) ));



