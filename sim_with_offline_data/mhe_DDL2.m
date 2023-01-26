function x_full = mhe_DDL2(x_old,beq,Aeq,H,f,n_states,n_meas,L)

x_old_hat = x_old(n_states+1:2*n_states,:);

% z_hat = quadprog(H.'*H,-x_old_hat.'*f.'*H,A,b,Aeq,beq);
z_hat = quadprog(H.'*H,-x_old_hat.'*f.'*H,[],[],Aeq,beq);

x_full = z_hat(end-n_meas*L-n_states*L+1:end-n_meas*L,:);