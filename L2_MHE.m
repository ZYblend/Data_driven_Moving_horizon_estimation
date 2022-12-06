function [x_hat,r] = L2_MHE(y,u,T,H0,H1,F,A)

x_hat0 = quadprog(2*(H0.'*H0),-2*H0.'*(y-H1*u));
% fun = @(x)norm(y-H0*x-H1*u);
% x_hat0 = fmincon(fun,zeros(10,1),[],[]);

x_hat = mpower(A,T)*x_hat0 + F*u;
r = y - H0*x_hat0 - H1*u;

