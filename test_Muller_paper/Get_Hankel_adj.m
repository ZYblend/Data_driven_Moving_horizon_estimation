function [H] = Get_Hankel_adj(alph,L,n)
% Calculate the multiplicative adjoint of the Hankel matrix
% Inputs:
%  - alph [N-by-1]: input vector 
%  - L [integer scalar]: Depth of the required Hankel matrix
%  - n [integer scalar]: spatial dimension of the ambient signal space
%
% Outputs:
% - H [m*L-by-N-L+1]: Resulting adjoint Hankel matrix of depth L


% Olugbenga Moses Anubi
% 1/11/2023


N_alph = length(alph);
H_alph = zeros(L,N_alph+L-1);

for col = 1:L
    H_alph(col,col+(0:N_alph-1)) = alph(:);
end

H = kron(H_alph,eye(n));


