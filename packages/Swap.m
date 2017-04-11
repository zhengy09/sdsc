function [ MC,Adj ] = Swap(MC,Adj,ni,nj)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    
    Temp = MC(:,ni);
    MC(:,ni) = MC(:,nj);
    MC(:,nj) = Temp;
    
    Temp2 = eye(size(Adj));
    Temp2(ni,ni) = 0;Temp2(nj,nj) = 0;
    Temp2(ni,nj) = 1;Temp2(nj,ni) = 1;
    
    Adj = Temp2*Adj*Temp2;
end

