 function x_full = FIE_L1(x0_old,ud,xd,yd,P,R,rho,L,n_states,n_meas,t_sample,u,y)


%% Unpacking inputs
% ud = zd{1,1};
% xd = zd{1,2};
% yd = zd{1,3};

%% Get parameters
[H_L1,H_L2,f2,Aeq2,~,~] = Get_L1MHE_Param(P,R,rho,t_sample,n_states,n_meas,ud,xd,yd);

%% call quadratic programming solver
beq = [u;y;zeros(n_states*t_sample,1)];

fun = @(z)(1/2)*norm(H_L2*z-f2*x0_old,2)^2 + (1/2)*norm(H_L1*z,1)^2;
z_hat = fmincon(fun,zeros(size(H_L2,2),1),[],[],Aeq2,beq);

x_full_t = z_hat(end-n_meas*t_sample-n_states*t_sample+1:end-n_meas*t_sample,:);

x_full = zeros(L*n_states,1);
x_full(1:n_states*t_sample) = x_full_t;



