function  [x0_hat,x_next] = mhe_MDL1(u,x_old,beq,Aeq,H_L1,H_L2,f,F,A,n_states,L)


fun = @(z)(1/2)*norm(H_L2*z-f*x_old,2)^2 + (1/2)*norm(H_L1*z,1)^2;
z_hat = fmincon(fun,zeros(size(H_L2,2),1),[],[],Aeq,beq);

x0_hat = z_hat(1:n_states);
x_next = mpower(A,L)*x0_hat + F*u;