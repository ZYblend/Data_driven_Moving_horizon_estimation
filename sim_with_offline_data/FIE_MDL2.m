 function [x0_hat, x_next] = FIE_MDL2(x0_old,A,B,C,P,R,rho,n_states,n_meas,t_sample,u,y)


%% Get parameters
[H,f,Aeq,H1,F] = Get_MDL2_MHE_Param(P,R,rho,t_sample,n_states,n_meas,A,B,C);

%% call quadratic programming solver
beq = y-H1*u;
z_hat = quadprog(H.'*H,-x0_old.'*f.'*H,[],[],Aeq,beq);
x0_hat = z_hat(1:n_states,:);
x_next = mpower(A,t_sample)*x0_hat + F*u;




