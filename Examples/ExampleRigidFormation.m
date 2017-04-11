
%% Example 1, rigid formation --> Franscico -- Identical decoupled 

clc;clear;close all;

N = 4;       % network size N*N
n = 2; m = 1; % dimension of state & input

%% generate communication graph
Gc = zeros(N^2);      %% communication graph
Dis = zeros(N^2,2);   %% positon of each node

for i = 1:N^2
    Dis(i,1) = mod(i,N); Dis(i,2) = (i-mod(i,N))/N+1;
    if mod(i,N) == 0
        Dis(i,1) = N;
        Dis(i,2) = (i-mod(i,N))/N;
    end
    if i > 0 && i < N                        %% bottom line
        Gc(i,i+1) = 1; Gc(i+1,i) = 1;
        Gc(i,i+N) = 1; Gc(i+N,i) = 1;
    elseif i > N*(N-1)+1 && i < N^2+1        %% upper line
        Gc(i-1,i) = 1; Gc(i,i-1) = 1;
    elseif mod(i,N) == 1 && i ~= 1           %% Left line
        Gc(i,i-N) = 1;   Gc(i-N,i) = 1;
    elseif mod(i,N) == 0 && i ~= N^2         %% Right line
        Gc(i,i+N) = 1;  Gc(N+i,i) = 1;
    else
        Gc(i,i+N) = 1; Gc(N+i,i) = 1;
        Gc(i,i-N) = 1; Gc(i-N,i) = 1;
        Gc(i,i+1) = 1; Gc(1+i,i) = 1;
        Gc(i,i-1) = 1; Gc(i-1,i) = 1;
    end       
end

Gp = zeros(N^2);    %% dynamic is seperate, no coupling
% Gex = chordalExt(Gc+eye(N^2));  

%% chordal extension
Gex = Gc;
for i = 1:N^2
    if i < N*(N-1) && mod(i,N) ~=0                        
        Gex(i,i + N+1) = 1;Gex(i+N+1,i) = 1;
    end       
end

Gex1 = Gex;
for i = 1:N^2
    if i < N*(N-1)+1 && mod(i,N) ~=1                        
        Gex1(i,i + N-1) = 1;Gex1(i+N-1,i) = 1;
    end       
end

Gex2 = Gc;
for i = 1:N^2
    if mod((i-mod(i,N))/N,2) == 0 && mod(i,N) ~= 0  && i < N*(N-1) +1                     
        Gex2(i,i + N+1) = 1;Gex2(i+N+1,i) = 1;
    end 
    if mod((i-mod(i,N))/N,2) == 1 && mod(i,N) ~= 1 && mod(i,N) ~= 0 && i < N*(N-1) +1                      
        Gex2(i,i + N-1) = 1;Gex2(i+N-1,i) = 1;
    end
    if mod(i,N) == 0 && mod((i-mod(i,N))/N,2) == 0 && i < N*(N-1) +1 
        Gex2(i,i + N-1) = 1;Gex2(i+N-1,i) = 1;
    end
end


[y,alpha] = checkIfChordal(Gex+eye(N^2))
[y,alpha] = checkIfChordal(Gex1+eye(N^2))
[y,alpha] = checkIfChordal(Gex2+eye(N^2))

Gex11 = chordalExt(Gex+eye(N^2));
Gex12 = chordalExt(Gex1+eye(N^2));
Gex13 = chordalExt(Gex1+eye(N^2));

%% Figure
ExampleGraphPlot(Gc,Dis,26)
ExampleGraphPlot(Gex,Dis)
ExampleGraphPlot(Gex11,Dis)
ExampleGraphPlot(Gex1,Dis)
ExampleGraphPlot(Gex12,Dis)
ExampleGraphPlot(Gex2,Dis)
ExampleGraphPlot(Gex13,Dis)

