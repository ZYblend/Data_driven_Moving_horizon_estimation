 function x_full = FIE_DDL2(x0_old,ud,xd,yd,P,R,rho,L,n_states,n_meas,t_sample,u,y)


%% Unpacking inputs
% ud = zd{1,1};
% xd = zd{1,2};
% yd = zd{1,3};

%% Get parameters
[H,f,Aeq,A_ineq,b_ineq] = Get_DDL2_MHE_Param(P,R,rho,t_sample,n_states,n_meas,ud,xd,yd);

%% call quadratic programming solver
beq = [u;y;zeros(n_states*t_sample,1)];
z_hat = quadprog(H.'*H,-x0_old.'*f.'*H,[],[],Aeq,beq);
x_full_t = z_hat(end-n_meas*t_sample-n_states*t_sample+1:end-n_meas*t_sample,:);

x_full = zeros(L*n_states,1);
x_full(1:n_states*t_sample) = x_full_t;



