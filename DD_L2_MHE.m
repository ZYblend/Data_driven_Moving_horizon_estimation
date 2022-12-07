% function [x_hat] = DD_L2_MHE(y,u,T,Aeq,n_states,E_MH,CE_MH,DE_MH,n_alpha)
function [x_hat] = DD_L2_MHE(y,u,T,H_pinv,n_states,E_MH,DE_MH)

y_tilde = [u;y;E_MH*y-DE_MH*u];

%% least-square 
z = H_pinv*y_tilde;
x_hat = z((T-1)*n_states+1:T*n_states);  
