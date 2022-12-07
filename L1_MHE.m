function [x_hat] = L1_MHE(y,u,T,H0,H1,F,A,n_states)

n_meas_T = size(y,1);

c = [zeros(n_states,1);  % coefficient of x
     ones(n_meas_T,1)   % coefficient of eta
     ];
A_in = [-H0, -eye(n_meas_T);
        H0, -eye(n_meas_T)];
b_in = [-y+H1*u;
         y-H1*u];

z = linprog(c,A_in,b_in);
x_hat0 = z(1:n_states);

x_hat = mpower(A,T)*x_hat0 + F*u;