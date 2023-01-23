function [H] = Get_Hankel(u,L)
% Calculate Hanker matrix
%
% Inputs:
%        u: [N-by-m] data trajectory (N data, data dimention is m)
%
% Yu Zheng, FSU
% 12/05/2022

[N,m] = size(u);
H = zeros(m*L,N-L+1);

for col = 1:N-L+1
    u_col = u(col:col+L-1,:).';
    H(:,col) = u_col(:);
end

%% check if it is persistently exciting
if rank(H) == size(H,1)
    disp("persistently exciting!");
else
    disp("Not persistenly exciting...");
    disp("The rank of Hankel matrix is:");
    disp(num2str(rank(H)));
    disp("The size of Hankel matrix is:");
    disp(num2str(size(H)));
end


