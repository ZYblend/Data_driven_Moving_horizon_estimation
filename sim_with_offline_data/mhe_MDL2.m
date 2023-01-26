function [x0_hat,x_next] = mhe_MDL2(u,x_old,beq,Aeq,H,f,n_states,F,A,L)


z_hat = quadprog(H.'*H,-x_old.'*f.'*H,[],[],Aeq,beq);

x0_hat = z_hat(1:n_states,:);
x_next = mpower(A,L)*x0_hat + F*u;

