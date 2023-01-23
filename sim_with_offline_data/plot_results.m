%% Plot results

%% Extracting values
time_vec  = out.logsout.getElement('x').Values.Time;

% State vectors
x         = out.logsout.getElement('x').Values.Data;
x_hat_L2  = out.logsout.getElement('x_hat_L2').Values.Data; %% DD_L2
x_hat_L1  = out.logsout.getElement('x_hat_L1').Values.Data; %% DD_L1
y  = out.logsout.getElement('y').Values.Data;
ya  = out.logsout.getElement('ya').Values.Data;
u = out.logsout.getElement('u').Values.Data;

save('Data1.mat','x','x_hat_L1','x_hat_L2','u','y','ya','T_start_attack','I','-v7.3');

%% Plotting
LW = 1.5;  % linewidth
FS = 15;   % font size

figure
subplot(2,2,3)
plot(time_vec(1:200),x(1:200,1),'LineWidth',2*LW)
hold on, plot(time_vec(1:200),x_hat_L1(1:200,1),'LineWidth',LW);
title('x_1')
ylabel('DD L1-L2 MHE')
xline(T_start_attack,'LineWidth',0.5*LW)
set(gca,'fontweight','bold','fontsize',12) 
set(gca,'LineWidth',LW)
legend('Actual','Estimated')

subplot(2,2,4)
plot(time_vec(1:200),x(1:200,2),'LineWidth',2*LW)
hold on, plot(time_vec(1:200),x_hat_L1(1:200,2),'LineWidth',LW);
title('x_2')
xline(T_start_attack,'LineWidth',0.5*LW)
set(gca,'fontweight','bold','fontsize',12) 
set(gca,'LineWidth',LW)
legend('Actual','Estimated')

subplot(2,2,1)
plot(time_vec(1:200),x(1:200,1),'LineWidth',2*LW)
hold on, plot(time_vec(1:200),x_hat_L2(1:200,1),'LineWidth',LW)
ylabel('DD L2 MHE')
title('x_1')
xline(T_start_attack,'LineWidth',0.5*LW)
set(gca,'fontweight','bold','fontsize',12) 
set(gca,'LineWidth',LW)
legend('Actual','Estimated')

subplot(2,2,2)
plot(time_vec(1:200),x(1:200,2),'LineWidth',2*LW)
hold on, plot(time_vec(1:200),x_hat_L2(1:200,2),'LineWidth',LW)
title('x_2')
xline(T_start_attack,'LineWidth',0.5*LW)
set(gca,'fontweight','bold','fontsize',12) 
set(gca,'LineWidth',LW)
legend('Actual','Estimated')
