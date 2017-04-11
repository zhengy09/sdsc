function [Q,R] = Vector2Matrix(y,u,v,u1,v1,N,n,m,Nc,Nc1)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    Q = zeros(N*n);
    R = zeros(N*m,N*n);
    for i = 1:Nc
        if i < Nc1+1
            Q(u(i),v(i)) = y(i); Q(v(i),u(i)) = y(i);
        else
            R(u1(i-Nc1),v1(i-Nc1)) = y(i); 
        end
    end
    
end

