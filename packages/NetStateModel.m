function [ Awhole, Bwhole] = NetStateModel(A,B,Gp)
% Convert the data into state space form
%   Detailed explanation goes here

    [N,temp] = size(Gp);     % Number of nodes in the graph
    [n,m] = size(B{1});      % dimensions of input and state in each node 

    Awhole = zeros(N*n); Bwhole = zeros(N*n,N*m);
    for i = 1:N
        for j = 1:N
            if i == j
                Awhole((i-1)*n+1:i*n,(j-1)*n+1:j*n) = A{i,j};
                Bwhole((i-1)*n+1:i*n,(j-1)*m+1:j*m) = B{i};
            elseif Gp(i,j) == 1
                Awhole((i-1)*n+1:i*n,(j-1)*n+1:j*n) = A{i,j};
            end
        end
    end
end

