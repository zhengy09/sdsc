function [Q1,R1,R2] = Covert2Matrix(y)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

I = ones(2);
Q1 = zeros(2);
R1 = zeros(2);
R2 = zeros(2);

% correspond to the Q = [Q1 0; 0 Q2];

[u,v] = find(triu(I));
for i = 1:3 % corresponging to the diagnol elements of Q
    Q1(u(i),v(i)) = y(i); Q1(v(i),u(i)) = y(i);
end

% correspond to the R = [R11 R12; R21 R22];

[u,v] = find(I);
for i = 4:7% corresponging to the diagnol elements of R
    R1(u(i-3),v(i-3)) = y(i); 
end

[u,v] = find(I);
for i = 8:11 % corresponging to the diagnol elements of Z
    R2(u(i-7),v(i-7)) = y(i); 
end

end

