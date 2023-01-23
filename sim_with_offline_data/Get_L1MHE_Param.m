function [H_L1,H_L2,f,Aeq,A_ineq,b_ineq] = Get_L1MHE_Param(P,R,rho,L,n_states,n_meas,ud,xd,yd)
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

%% Get hankel matrix of order L
HL_u = Get_Hankel(ud,L);
HL_x = Get_Hankel(xd,L);
HL_y = Get_Hankel(yd,L);


%% Get objective parameters
n_alpha = size(HL_u,2);

H1 = zeros(1,n_alpha+n_states*L + n_meas*L);

H21 = [zeros(n_states,n_alpha) (sqrt(rho)^L)*sqrt(P)*eye(n_states) zeros(n_states,n_states*(L-1)) zeros(n_states,n_meas*L)];
H22 = zeros(n_states*(L-1),n_alpha+n_states*L+ n_meas*L);
H2 = [H21; H22];

H3 = zeros(n_meas*L,n_alpha+n_states*L++ n_meas*L);

H_L2 = [H1;
     H2;
     H3];

f = [zeros(1,n_states);
     eye(n_states);
     zeros(n_states*(L-1),n_states);
     zeros(n_meas*L,n_states)];

rho_forget = zeros(n_meas*L);
for idx = 1:L
    rho_forget(n_meas*(idx-1)+1:n_meas*idx,n_meas*(idx-1)+1:n_meas*idx) = sqrt(rho)^(L+1-idx)*eye(n_meas);
end

e_select = R*[zeros(n_meas*L,n_alpha) zeros(n_meas*L,n_states*L) eye(n_meas*L)];
H_L1 = rho_forget*e_select;

%% Get constraint parameters
Aeq = [HL_u zeros(size(HL_u,1),n_states*L) zeros(size(HL_u,1),n_meas*L);
       HL_y zeros(size(HL_y,1),n_states*L) eye(n_meas*L);
       HL_x -eye(n_states*L) zeros(n_states*L,n_meas*L)];

A_ineq = -[zeros(1,n_alpha) zeros(1,n_states*L) zeros(1,n_meas*L);
          zeros(n_states*L,n_alpha) eye(n_states*L) zeros(n_states*L,n_meas*L);
          zeros(n_meas*L,n_alpha) zeros(n_meas*L,n_states*L) zeros(n_meas*L)];
b_ineq = zeros(1+n_states*L+n_meas*L,1);