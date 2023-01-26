function  x_full = mhe_DDL1(x_old,beq,Aeq,H_L1,H_L2,f,n_states,n_meas,L)

x_old_hat = x_old(n_states+1:2*n_states,:);

fun = @(z)(1/2)*norm(H_L2*z-f*x_old_hat,2)^2 + (1/2)*norm(H_L1*z,1)^2;
z_hat = fmincon(fun,zeros(size(H_L2,2),1),[],[],Aeq,beq);

x_full = z_hat(end-n_meas*L-n_states*L+1:end-n_meas*L,:);