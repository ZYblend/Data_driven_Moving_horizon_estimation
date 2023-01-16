function [H] = Get_Hankel(u,L)
% Calculate Hanker matrix
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
if rank(H) == m*L
    disp("persistently exciting!");
else
    disp("Not persistenly exciting...");
    disp("The rank of hanker matrix is:");
    disp(num2str(rank(H)));
end


