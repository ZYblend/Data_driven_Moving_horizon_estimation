% function [x_hat] = DD_L2_MHE(y,u,T,Aeq,n_states,E_MH,CE_MH,DE_MH,n_alpha)
function [x_hat] = DD_L2_MHE(y,u,T,Aeq,n_states,E_MH,CE_MH,DE_MH,n_alpha)

% x = [x_hat; e; alpha]
n_meas_T = size(E_MH,2);

beq = [u;y;E_MH*y-DE_MH*u];

%% least-square (not work)
% xhat = Aeq_inv*beq;
% xhat = quadprog(2*(Aeq.'*Aeq),-2*Aeq.'*beq);

%% equality constraints (work but slow)
%% function [x_hat] = DD_L2_MHE(y,u,T,Aeq,n_states,E_MH,CE_MH,DE_MH,n_alpha)
H = diag([zeros(n_states*T,1); ones(n_meas_T,1); zeros(n_alpha,1)]);
xhat = quadprog((H.'*H),zeros(1,n_states*T+n_meas_T+n_alpha),[],[],Aeq,beq);


x_hat = xhat((T-1)*n_states+1:T*n_states);  
% 
% r = E_MH*y - CE_MH*xhat(1:T*n_states) - DE_MH*u;