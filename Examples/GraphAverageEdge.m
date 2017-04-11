function Gc = GraphAverageEdge(G,Gp);
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

p = 0.1;

Gtemp = (G+Gp)/2;
[u,v] = find(Gtemp == 1/2);

Num = length(u);

Index = randi(Num, floor(Num*p)+1,1);
    Gc = Gp;
    Gc(u(Index),v(Index)) = 1;



end

