function [x_hat] = DD_L1_MHE(y,u,T,A_in1, A_in2,n_states,E_MH,CE_MH,DE_MH,n_alpha)

% x = [x_hat; e; alpha]
n_meas_T = size(E_MH,2);

y_tilde = [u;y;E_MH*y-DE_MH*u];

%% equality constraints
%% function [x_hat] = DD_L2_MHE(y,u,T,Aeq,n_states,E_MH,CE_MH,DE_MH,n_alpha)
c = [zeros(n_states*T+n_alpha,1);  ones(n_meas_T,1)];
A_in = [-A_in1, -A_in2;
        A_in1, -A_in2];
b_in = [-y_tilde;
         y_tilde];
z = linprog(c,A_in,b_in);

x_hat = z((T-1)*n_states+1:T*n_states);  
% 
% r = E_MH*y - CE_MH*xhat(1:T*n_states) - DE_MH*u;