
function [ Adj1,MC1 ] = TreeRegulation( Adj,MC )
% Adj, MC

[n1,n1] = size(Adj);
Adj1 = Adj;
MC1 = MC;
for i = 1:n1
    Atemp = triu(Adj1);
    Index = find(Atemp(i,:) == 1);
    if i == 1 
        for j = 1:length(Index)
            [MC1,Adj1] = Swap(MC1,Adj1,j+1,Index(j));
        end
    end
    if i > 1
        [uu,vv] = find(Atemp(1:i-1,:) == 1);
        for j = 1:length(Index)
            [MC1,Adj1] = Swap(MC1,Adj1,vv(end)+j,Index(j));
        end
        
    end
end
end