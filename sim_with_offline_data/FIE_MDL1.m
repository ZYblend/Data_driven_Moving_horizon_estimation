 function [x0_hat, x_next] = FIE_MDL1(x0_old,A,B,C,P,R,rho,n_states,n_meas,t_sample,u,y)
%% Model-driven resilient MHE 

%% Get parameters
[H_L1,H_L2,f,Aeq,H1,F] = Get_MDL1_MHE_Param(P,R,rho,t_sample,n_states,n_meas,A,B,C);

%% call nonlinear programming solver
beq = y-H1*u;

fun = @(z)(1/2)*norm(H_L2*z-f*x0_old,2)^2 + (1/2)*norm(H_L1*z,1)^2;
z_hat = fmincon(fun,zeros(size(H_L2,2),1),[],[],Aeq,beq);

x0_hat = z_hat(1:n_states);

x_next = mpower(A,t_sample)*x0_hat + F*u;



