function [x_hat] = DD_L1_MHE(y,u,T,A_in1, A_in2,n_states,E_MH,DE_MH,n_alpha)


n_meas_T = size(E_MH,2);
y_tilde = [u;y;E_MH*y-DE_MH*u];

% %% linear programming
c = [zeros(n_states*T+n_alpha,1);  ones(n_meas_T,1)];
A_in = [-A_in1, -A_in2;
        A_in1, -A_in2];
b_in = [-y_tilde;
         y_tilde];
z = linprog(c,A_in,b_in);

x_hat = z((T-1)*n_states+1:T*n_states);  

% %% fmincon
% fun = @(x)c.'*x;
% Aeq = [A_in1,A_in2];
% z = fmincon(fun,zeros(n_states*T+n_alpha+n_meas_T,1),[],[],Aeq,y_tilde);
% x_hat = z((T-1)*n_states+1:T*n_states);  