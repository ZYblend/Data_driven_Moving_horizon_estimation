%% Plot error table
% x axis: number of attacks
% y axis: different MHE scheme: MD_L2, DD_l2, MD_L1, DD_L1
% Two error metrics will be calculated: root mean square error, maximum absolute error
%
% Yu Zheng, Florida State University
% 12/13/2022


metric_rms = @(x,dim) sqrt(sum(x.^2,dim)/size(x,dim));
metric_max = @(x,dim) max(abs(x),[],dim);

T_sample = 0.01;

%% number of attacks = 1
x = load('Results/Data1').x;
x_hat1 = load('Results/Data1').x_hat1;
x_hat2 = load('Results/Data1').x_hat2;
x_hat3 = load('Results/Data1').x_hat3;
x_hat4 = load('Results/Data1').x_hat4;
T_start_opt = load('Results/Data1').T_start_opt;
N_start_opt = round(T_start_opt/T_sample);

% calculate error metrices
e_MD_L2 = vecnorm(x(N_start_opt+1:end,:) - x_hat1(N_start_opt+1:end,:),2,2);
e_DD_L2 = vecnorm(x(N_start_opt+1:end,:) - x_hat2(N_start_opt+1:end,:),2,2);
e_MD_L1 = vecnorm(x(N_start_opt+1:end,:) - x_hat3(N_start_opt+1:end,:),2,2);
e_DD_L1 = vecnorm(x(N_start_opt+1:end,:) - x_hat4(N_start_opt+1:end,:),2,2);

rms_MD_L2_1 = metric_rms(e_MD_L2,1);
rms_DD_L2_1 = metric_rms(e_DD_L2,1);
rms_MD_L1_1 = metric_rms(e_MD_L1,1);
rms_DD_L1_1 = metric_rms(e_DD_L1,1);

max_MD_L2_1 = metric_max(e_MD_L2,1);
max_DD_L2_1 = metric_max(e_DD_L2,1);
max_MD_L1_1 = metric_max(e_MD_L1,1);
max_DD_L1_1 = metric_max(e_DD_L1,1);


%% number of attacks = 2
x = load('Results/Data2.mat').x;
x_hat1 = load('Results/Data2.mat').x_hat1;
x_hat2 = load('Results/Data2.mat').x_hat2;
x_hat3 = load('Results/Data2.mat').x_hat3;
x_hat4 = load('Results/Data2.mat').x_hat4;
T_start_opt = load('Results/Data2.mat').T_start_opt;
N_start_opt = round(T_start_opt/T_sample);

% calculate error metrices
e_MD_L2 = vecnorm(x(N_start_opt+1:end,:) - x_hat1(N_start_opt+1:end,:),2,2);
e_DD_L2 = vecnorm(x(N_start_opt+1:end,:) - x_hat2(N_start_opt+1:end,:),2,2);
e_MD_L1 = vecnorm(x(N_start_opt+1:end,:) - x_hat3(N_start_opt+1:end,:),2,2);
e_DD_L1 = vecnorm(x(N_start_opt+1:end,:) - x_hat4(N_start_opt+1:end,:),2,2);

rms_MD_L2_2 = metric_rms(e_MD_L2,1);
rms_DD_L2_2 = metric_rms(e_DD_L2,1);
rms_MD_L1_2 = metric_rms(e_MD_L1,1);
rms_DD_L1_2 = metric_rms(e_DD_L1,1);

max_MD_L2_2 = metric_max(e_MD_L2,1);
max_DD_L2_2 = metric_max(e_DD_L2,1);
max_MD_L1_2 = metric_max(e_MD_L1,1);
max_DD_L1_2 = metric_max(e_DD_L1,1);

%% number of attacks = 2
x = load('Results/Data3.mat').x;
x_hat1 = load('Results/Data3.mat').x_hat1;
x_hat2 = load('Results/Data3.mat').x_hat2;
x_hat3 = load('Results/Data3.mat').x_hat3;
x_hat4 = load('Results/Data3.mat').x_hat4;
T_start_opt = load('Results/Data3.mat').T_start_opt;
N_start_opt = round(T_start_opt/T_sample);

% calculate error metrices
e_MD_L2 = vecnorm(x(N_start_opt+1:end,:) - x_hat1(N_start_opt+1:end,:),2,2);
e_DD_L2 = vecnorm(x(N_start_opt+1:end,:) - x_hat2(N_start_opt+1:end,:),2,2);
e_MD_L1 = vecnorm(x(N_start_opt+1:end,:) - x_hat3(N_start_opt+1:end,:),2,2);
e_DD_L1 = vecnorm(x(N_start_opt+1:end,:) - x_hat4(N_start_opt+1:end,:),2,2);

rms_MD_L2_3 = metric_rms(e_MD_L2,1);
rms_DD_L2_3 = metric_rms(e_DD_L2,1);
rms_MD_L1_3 = metric_rms(e_MD_L1,1);
rms_DD_L1_3 = metric_rms(e_DD_L1,1);

max_MD_L2_3 = metric_max(e_MD_L2,1);
max_DD_L2_3 = metric_max(e_DD_L2,1);
max_MD_L1_3 = metric_max(e_MD_L1,1);
max_DD_L1_3 = metric_max(e_DD_L1,1);


%% plot table
rms_table = [rms_MD_L2_1, rms_MD_L2_2, rms_MD_L2_3;
             rms_DD_L2_1, rms_DD_L2_2, rms_DD_L2_3;
             rms_MD_L1_1, rms_MD_L1_2, rms_MD_L1_3;
             rms_DD_L1_1, rms_DD_L1_2, rms_DD_L1_3]

max_table = [max_MD_L2_1, max_MD_L2_2, max_MD_L2_3;
             max_DD_L2_1, max_DD_L2_2, max_DD_L2_3;
             max_MD_L1_1, max_MD_L1_2, max_MD_L1_3;
             max_DD_L1_1, max_DD_L1_2, max_DD_L1_3]
