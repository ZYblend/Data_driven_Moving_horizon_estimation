function [H] = Get_Hanker(u,L)

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
    disp(num2str(rank(H)));
end


