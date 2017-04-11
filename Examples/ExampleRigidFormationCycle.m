
%% Example 1, rigid formation --> Franscico -- Identical decoupled 

clc;clear;close all;

N = 16;       % network size N*N
n = 2; m = 1; % dimension of state & input

%% generate communication graph && banded communication
Gc = zeros(N);      %% communication graph
Dis = zeros(N,2);   %% positon of each node
Bandwith = 2;

Angle = -pi/2+pi/N:2*pi/N:pi/2*3-pi/N; 
R = 100;
for i = 1:N
    Dis(i,1) = R*cos(Angle(i));Dis(i,2) = R*sin(Angle(i));
    if i < Bandwith+1                        %% bottom line
        for j = 1:i+Bandwith
            Gc(i,j) = 1; Gc(j,i) = 1;
        end
    elseif i > N - Bandwith-1
        for j = i-Bandwith:N
            Gc(i,j) = 1; Gc(j,i) = 1;
        end
    else
        for j = i-Bandwith:i+Bandwith
            Gc(i,j) = 1; Gc(j,i) = 1;
        end
    end       
end

Gp = zeros(N^2);    %% dynamic is seperate, no coupling
% Gex = chordalExt(Gc+eye(N^2));  



[y,alpha] = checkIfChordal(Gc)

%% Figure
ExampleGraphPlot(Gc,Dis,26)
ExampleGraphPlot(Gp,Dis,26)
Kc = kron(Gc+eye(size(Gc)),ones(1,2));
figure
set(gcf,'position',[300,200,550,340])
spy(Kc)



