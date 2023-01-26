function [H,f,Aeq,H1,F] = Get_MDL2_MHE_Param(P,R,rho,L,n_states,n_meas,A,B,C)
%% This function is to get the vectorized parameter for MHE and FIE program
% Inputs:
%        - P,R: weights for terminal cost and stage cost respectively
%        - rho: forget factor
%        - L: horizon length
%        - n_states, n_meas: dimention of system
%        - ud, xd, yd: historical data (offline)
%
% Output: 
%        Parameters for MHE: z = [alpha^T, xhat^T, sigma^T]^T
%                Minimize (1/2)||Hz-f*xhat_old||_2
%               Subject to Aeq z = beq
%                          A_ineq z <= b_ineq
%
% Yu Zheng. Florida State University
% 01/18/2023

[H0,H1,F] = opti_params(A,B,C,L);
% H0: state-ouput linear map                       [n_meas*T-by-n_states]
% H1:  input-output linear map                     [n_meas*T-by-n_int*(T-1)]
% F:  Observer input-state propagation matrix      [n_meas-by-n_int*(T-1)]

%% Get objective parameters

H2 = [(sqrt(rho)^L)*sqrt(P)*eye(n_states) zeros(n_states,n_meas*L)];

rho_forget = zeros(n_meas*L);
for idx = 1:L
    rho_forget(n_meas*(idx-1)+1:n_meas*idx,n_meas*(idx-1)+1:n_meas*idx) = sqrt(rho)^(L+1-idx)*eye(n_meas);
end
H3 = [zeros(n_meas*L,n_states) rho_forget*sqrt(R)];

H = [H2;
     H3];

f = [eye(n_states);
     zeros(n_meas*L,n_states)];

%% Get constraint parameters
Aeq = [H0 eye(n_meas*L)];
