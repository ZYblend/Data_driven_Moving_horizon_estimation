%% Plot error table
% x axis: number of attacks
% y axis: different MHE scheme: MD_L2, DD_l2, MD_L1, DD_L1
% Two error metrics will be calculated: root mean square error, maximum absolute error
%
% Yu Zheng, Florida State University
% 12/13/2022

%% run all sims
% comment line 115, 116 in Run_model.m: 115 n_attack =  round(0.2*n_meas); 
%                                       116 I = randperm(n_meas,n_attack);
% 1. no attack
I = [];
Run_model;
out = sim("test_sys_MHE.slx");
x           = out.logsout.getElement('x').Values.Data;
x_hat_DDL2  = out.logsout.getElement('x_hat_DDL2').Values.Data; %% DD_L2
x_hat_DDL1  = out.logsout.getElement('x_hat_DDL1').Values.Data; %% DD_L1
x_hat_MDL2  = out.logsout.getElement('x_hat_MDL2').Values.Data; %% MD_L2
x_hat_MDL1  = out.logsout.getElement('x_hat_MDL1').Values.Data; %% MD_L1

filename = "results/Data_0_attacks_1.mat";

save(filename,'x','x_hat_DDL1','x_hat_DDL2','x_hat_MDL1','x_hat_MDL2','T_start_attack','I','-v7.3');

% 2. 1 attack
for iter =1:6
    I = iter;
    Run_model;
    out = sim("test_sys_MHE.slx");
    x           = out.logsout.getElement('x').Values.Data;
    x_hat_DDL2  = out.logsout.getElement('x_hat_DDL2').Values.Data; %% DD_L2
    x_hat_DDL1  = out.logsout.getElement('x_hat_DDL1').Values.Data; %% DD_L1
    x_hat_MDL2  = out.logsout.getElement('x_hat_MDL2').Values.Data; %% MD_L2
    x_hat_MDL1  = out.logsout.getElement('x_hat_MDL1').Values.Data; %% MD_L1

    filename = "results/Data_1_attacks_"+num2str(iter)+".mat";

    save(filename,'x','x_hat_DDL1','x_hat_DDL2','x_hat_MDL1','x_hat_MDL2','T_start_attack','I','-v7.3');
end

% 2. 2 attacks
I_2_attacks = nchoosek(1:6,2);
for iter =1:15
    I = I_2_attacks(iter,:);
    Run_model;
    out = sim("test_sys_MHE.slx");
    x           = out.logsout.getElement('x').Values.Data;
    x_hat_DDL2  = out.logsout.getElement('x_hat_DDL2').Values.Data; %% DD_L2
    x_hat_DDL1  = out.logsout.getElement('x_hat_DDL1').Values.Data; %% DD_L1
    x_hat_MDL2  = out.logsout.getElement('x_hat_MDL2').Values.Data; %% MD_L2
    x_hat_MDL1  = out.logsout.getElement('x_hat_MDL1').Values.Data; %% MD_L1

    filename = "results/Data_2_attacks_"+num2str(iter)+".mat";

    save(filename,'x','x_hat_DDL1','x_hat_DDL2','x_hat_MDL1','x_hat_MDL2','T_start_attack','I','-v7.3');
end

% 2. 3 attacks
I_3_attacks = nchoosek(1:6,3);
for iter =1:20
    I = I_3_attacks(iter,:);
    Run_model;
    out = sim("test_sys_MHE.slx");
    x           = out.logsout.getElement('x').Values.Data;
    x_hat_DDL2  = out.logsout.getElement('x_hat_DDL2').Values.Data; %% DD_L2
    x_hat_DDL1  = out.logsout.getElement('x_hat_DDL1').Values.Data; %% DD_L1
    x_hat_MDL2  = out.logsout.getElement('x_hat_MDL2').Values.Data; %% MD_L2
    x_hat_MDL1  = out.logsout.getElement('x_hat_MDL1').Values.Data; %% MD_L1

    filename = "results/Data_3_attacks_"+num2str(iter)+".mat";

    save(filename,'x','x_hat_DDL1','x_hat_DDL2','x_hat_MDL1','x_hat_MDL2','T_start_attack','I','-v7.3');
end


%% plot table
metric_rms = @(x,dim) sqrt(sum(x.^2,dim)/size(x,dim));
metric_max = @(x,dim) max(abs(x),[],dim);

T_sample = 0.01;

tot_attacks = 3;
rms_metric = zeros(4,tot_attacks+1);
max_metric = zeros(4,tot_attacks+1);
num_cases = [1,6,15,20];

for n_attack = 0:tot_attacks
    for idx = 1:num_cases(n_attack+1)
        filename = "results/Data_"+num2str(n_attack)+"_attacks_"+num2str(idx)+".mat";
        x = load(filename).x;
        x_hat_DDL2 = load(filename).x_hat_DDL2;
        x_hat_DDL1 = load(filename).x_hat_DDL1;
        x_hat_MDL2 = load(filename).x_hat_MDL2;
        x_hat_MDL1 = load(filename).x_hat_MDL1;
    
    
        % calculate error metrices
        e_DDL2 =  vecnorm(x - x_hat_DDL2,2,2);
        e_DDL1 = vecnorm(x - x_hat_DDL1,2,2);
        e_MDL2 = vecnorm(x - x_hat_MDL2,2,2);
        e_MDL1 = vecnorm(x - x_hat_MDL1,2,2);
        
        rms_metric(2,n_attack+1) = rms_metric(2,n_attack+1) + metric_rms(e_DDL2,1);
        rms_metric(4,n_attack+1) = rms_metric(4,n_attack+1) + metric_rms(e_DDL1,1);
        rms_metric(1,n_attack+1) = rms_metric(1,n_attack+1) + metric_rms(e_MDL2,1);
        rms_metric(3,n_attack+1) = rms_metric(3,n_attack+1) + metric_rms(e_MDL1,1);
        
        max_metric(2,n_attack+1) = max(max_metric(2,n_attack+1), metric_max(e_DDL2,1));
        max_metric(4,n_attack+1) = max(max_metric(4,n_attack+1), metric_max(e_DDL1,1));
        max_metric(1,n_attack+1) = max(max_metric(1,n_attack+1), metric_max(e_MDL2,1));
        max_metric(3,n_attack+1) = max(max_metric(3,n_attack+1), metric_max(e_MDL1,1));
    end
    rms_metric(:,n_attack+1) = rms_metric(:,n_attack+1)/num_cases(n_attack+1);
end

%% plot table
Approach = {'Robust MD-MHE'; 'Robust DD-MHE';'Resilient MD-MHE'; 'Resilient DD-MHE'};
zero_attack = rms_metric(:,1);
one_attack = rms_metric(:,2);
two_attacks = rms_metric(:,3);
three_attacks = rms_metric(:,4);
rms_table = table(Approach, zero_attack, one_attack, two_attacks, three_attacks)

zero_attack = max_metric(:,1);
one_attack = max_metric(:,2);
two_attacks = max_metric(:,3);
three_attacks = max_metric(:,4);
max_table = table(Approach, zero_attack, one_attack, two_attacks, three_attacks)
