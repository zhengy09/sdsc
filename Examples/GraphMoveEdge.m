function [ Gn ] = GraphMoveEdge(G, p)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    [u,v] = find(G == 1);
    Num = length(u);
    Index = randi(Num, floor(Num*p)+1,1);
    Gn = G;
    Gn(u(Index),v(Index)) = 0;


end

