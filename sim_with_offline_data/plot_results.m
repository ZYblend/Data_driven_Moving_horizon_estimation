%% Plot results

%% Extracting values
time_vec  = out.logsout.getElement('x').Values.Time;

% State vectors
x         = out.logsout.getElement('x').Values.Data;
x_hat_DDL2  = out.logsout.getElement('x_hat_DDL2').Values.Data; %% DD_L2
x_hat_DDL1  = out.logsout.getElement('x_hat_DDL1').Values.Data; %% DD_L1
x_hat_MDL2  = out.logsout.getElement('x_hat_MDL2').Values.Data; %% MD_L2
x_hat_MDL1  = out.logsout.getElement('x_hat_MDL1').Values.Data; %% MD_L1
y  = out.logsout.getElement('y').Values.Data;
ya  = out.logsout.getElement('ya').Values.Data;
u = out.logsout.getElement('u').Values.Data;

save('results/Data_3_attacks.mat','x','x_hat_DDL1','x_hat_DDL2','x_hat_MDL1','x_hat_MDL2','u','y','ya','T_start_attack','I','-v7.3');

%% Plotting
LW = 1.5;  % linewidth
FS = 25;   % font size

figure
% x1
% MD L2
subplot(4,2,1)
plot(time_vec(1:200),x(1:200,1),'LineWidth',2*LW)
hold on, plot(time_vec(1:200),x_hat_MDL2(1:200,1),'LineWidth',LW);
title('x_1')
xline(T_start_attack,'LineWidth',0.5*LW)
set(gca,'fontweight','bold','fontsize',FS) 
set(gca,'LineWidth',LW)
ylabel('Robust MD-MHE','fontsize',22)
% legend('Actual','Estimated')

% DD L2
subplot(4,2,3)
plot(time_vec(1:200),x(1:200,1),'LineWidth',2*LW)
hold on, plot(time_vec(1:200),x_hat_DDL2(1:200,1),'LineWidth',LW);
xline(T_start_attack,'LineWidth',0.5*LW)
set(gca,'fontweight','bold','fontsize',FS) 
set(gca,'LineWidth',LW)
ylabel('Robust DD-MHE','fontsize',22)
% legend('Actual','Estimated')

% MD L1
subplot(4,2,5)
plot(time_vec(1:200),x(1:200,1),'LineWidth',2*LW)
hold on, plot(time_vec(1:200),x_hat_MDL1(1:200,1),'LineWidth',LW)
xline(T_start_attack,'LineWidth',0.5*LW)
set(gca,'fontweight','bold','fontsize',FS) 
set(gca,'LineWidth',LW)
ylabel('Resilient MD-MHE','fontsize',22)
% legend('Actual','Estimated')

% DD L1
subplot(4,2,7)
plot(time_vec(1:200),x(1:200,1),'LineWidth',2*LW)
hold on, plot(time_vec(1:200),x_hat_DDL1(1:200,1),'LineWidth',LW)
xline(T_start_attack,'LineWidth',0.5*LW)
set(gca,'fontweight','bold','fontsize',FS) 
set(gca,'LineWidth',LW)
ylabel('Resilient DD-MHE','fontsize',22')
% legend('Actual','Estimated')

% x2
% MD L2
subplot(4,2,2)
plot(time_vec(1:200),x(1:200,2),'LineWidth',2*LW)
hold on, plot(time_vec(1:200),x_hat_MDL2(1:200,2),'LineWidth',LW);
title('x_2')
% ylabel('DD L1-L2 MHE')
xline(T_start_attack,'LineWidth',0.5*LW)
set(gca,'fontweight','bold','fontsize',FS) 
set(gca,'LineWidth',LW)
% legend('Actual','Estimated')

% DD L2
subplot(4,2,4)
plot(time_vec(1:200),x(1:200,2),'LineWidth',2*LW)
hold on, plot(time_vec(1:200),x_hat_DDL2(1:200,2),'LineWidth',LW);
xline(T_start_attack,'LineWidth',0.5*LW)
set(gca,'fontweight','bold','fontsize',FS) 
set(gca,'LineWidth',LW)
% legend('Actual','Estimated')

% MD L1
subplot(4,2,6)
plot(time_vec(1:200),x(1:200,2),'LineWidth',2*LW)
hold on, plot(time_vec(1:200),x_hat_MDL1(1:200,2),'LineWidth',LW)
xline(T_start_attack,'LineWidth',0.5*LW)
set(gca,'fontweight','bold','fontsize',FS) 
set(gca,'LineWidth',LW)
% legend('Actual','Estimated')

% DD L1
subplot(4,2,8)
plot(time_vec(1:200),x(1:200,2),'LineWidth',2*LW)
hold on, plot(time_vec(1:200),x_hat_DDL1(1:200,2),'LineWidth',LW)
xline(T_start_attack,'LineWidth',0.5*LW)
set(gca,'fontweight','bold','fontsize',FS) 
set(gca,'LineWidth',LW)
% legend('Actual','Estimated')
