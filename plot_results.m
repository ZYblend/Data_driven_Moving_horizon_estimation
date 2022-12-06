%% Plot results

%% Extracting values
time_vec  = out.logsout.getElement('x').Values.Time;

% State vectors
x         = out.logsout.getElement('x').Values.Data;
x_hat  = out.logsout.getElement('x_hat').Values.Data; 
y  = out.logsout.getElement('y').Values.Data;
ya  = out.logsout.getElement('ya').Values.Data;
u  = out.logsout.getElement('u').Values.Data;

% save('MD_L2.mat','x','x_hat','y','ya','u','T_start_opt', 'T_start_attack' ,'I','-v7.3');
% save('DD_L2.mat','x','x_hat','y','ya','u','T_start_opt', 'T_start_attack' ,'I','-v7.3');
% save('MD_L1.mat','x','x_hat','y','ya','u','T_start_opt', 'T_start_attack','I' ,'-v7.3');
save('DD_L1.mat','x','x_hat','y','ya','u','T_start_opt', 'T_start_attack' ,'-v7.3');

plot(x-x_hat)