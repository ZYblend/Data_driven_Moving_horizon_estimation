%% Plot results

%% Extracting values
time_vec  = out.logsout.getElement('x').Values.Time;

% State vectors
x         = out.logsout.getElement('x').Values.Data;
x_hat1  = out.logsout.getElement('x_hat1').Values.Data; %% MD_L2
x_hat2  = out.logsout.getElement('x_hat2').Values.Data; %% DD_L2
x_hat3  = out.logsout.getElement('x_hat3').Values.Data; %% MD_L1
x_hat4  = out.logsout.getElement('x_hat4').Values.Data; %% DD_L1
y  = out.logsout.getElement('y').Values.Data;
ya  = out.logsout.getElement('ya').Values.Data;


save('Data1.mat','x','x_hat1','x_hat2','x_hat3','x_hat4','y','ya','T_start_opt', 'T_start_attack','I','-v7.3');

%% Plotting
LW = 1.5;  % linewidth
FS = 15;   % font size

figure
subplot(4,2,1)
plot(x(:,1),'LineWidth',LW)
hold on, plot(x_hat1(:,1),'LineWidth',LW);
title('x_1')
ylabel('MD L2 MHE')
xlim([0,800])
ylim([-40,20])
set(gca,'fontweight','bold','fontsize',12) 
set(gca,'LineWidth',LW)

subplot(4,2,2)
plot(x(:,2),'LineWidth',LW)
hold on, plot(x_hat1(:,2),'LineWidth',LW);
title('x_2')
xlim([0,800])
ylim([-40,20])
set(gca,'fontweight','bold','fontsize',12) 
set(gca,'LineWidth',LW)

subplot(4,2,3)
plot(x(:,1),'LineWidth',LW)
hold on, plot(x_hat2(:,1),'LineWidth',LW)
ylabel('DD L2 MHE')
xlim([0,800])
ylim([-40,20])
set(gca,'fontweight','bold','fontsize',12) 
set(gca,'LineWidth',LW)

subplot(4,2,4)
plot(x(:,2),'LineWidth',LW)
hold on, plot(x_hat2(:,2),'LineWidth',LW)
xlim([0,800])
ylim([-40,20])
set(gca,'fontweight','bold','fontsize',12) 
set(gca,'LineWidth',LW)

subplot(4,2,5)
plot(x(:,1),'LineWidth',LW)
hold on, plot(x_hat3(:,1),'LineWidth',LW)
ylabel('MD L1 MHE')
xlim([0,800])
ylim([-40,20])
set(gca,'fontweight','bold','fontsize',12) 
set(gca,'LineWidth',LW)

subplot(4,2,6)
plot(x(:,2),'LineWidth',LW)
hold on, plot(x_hat3(:,2),'LineWidth',LW)
xlim([0,800])
ylim([-40,20])
set(gca,'fontweight','bold','fontsize',12) 
set(gca,'LineWidth',LW)

subplot(4,2,7)
plot(x(:,1),'LineWidth',LW)
hold on, plot(x_hat4(:,1),'LineWidth',LW)
ylabel('DD L1 MHE')
xlim([0,800])
ylim([-40,20])
set(gca,'fontweight','bold','fontsize',12) 
set(gca,'LineWidth',LW)


subplot(4,2,8)
plot(x(:,2),'LineWidth',LW)
hold on, plot(x_hat4(:,2),'LineWidth',LW)
xlim([0,800])
ylim([-40,20])
set(gca,'fontweight','bold','fontsize',12) 
set(gca,'LineWidth',LW)