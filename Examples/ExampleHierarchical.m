
%% Example 1.1 hierarchical system in Fig.3
clc;clear;close all
%% Data Generation Part
N = 8;     %% the number of nodes in the graph
n = 2;     %% the dimension of state in each node
m = 1;     %% the dimension of input in each node

%% Dynamics
Gp = [0 0 0 0 0 0 0 0;
      1 0 0 0 0 0 0 0;
      1 0 0 0 0 0 0 0;
      1 0 0 0 0 0 0 0;
      0 1 0 0 0 0 0 0;
      0 1 1 0 0 0 0 0;
      0 0 1 1 0 0 0 0;
      0 0 0 1 0 0 0 0;];
Gc = Gp';



% Dynamic matrices
A = cell(N,N);   % matrices for A
B = cell(N,1);   % matrices for B

% Dynamical part
for i = 1:N
    for j = 1:N
        if Gp(i,j) == 1
            A{i,j} = exp(-(i-j)^2/10)*eye(n);     %% node j has influence on node i
        end  
        if i == j
            B{i} = [0;1];
            A{i,i} = [1 1; 1 2];
        end      
    end
end

%  load dataFailure
[Awhole, Bwhole] = NetStateModel(A,B,Gp);


%% Design of structred feedback gains I: centralized solution (block diagonal Q)
[K,Q,R,SedumiInfo,Tdata,Tsolver] = scCentr(A,B,Gp,Gc);
[Tdata,Tsolver,Tdata+Tsolver]

Flag = StrucCheck(K,Gc,m,n);
if sum(real(eig(Awhole+Bwhole*K)) > 0) == 0
    disp('******** The closed-loop system is stable. ********')
end

figure;
set(gcf,'position',[300,200,350,200])
spy(K)


%% Design of structred feedback gains II: Sequential design method
%% Note that this will result undesesirable large gains
[Kse,Qse,Rse,Tdatase,Tsolverse,Tgraph,Ttotal] = scSeq(A,B,Gp,Gc);

Flag = StrucCheck(Kse,Gc,m,n);
if sum(real(eig(Awhole+Bwhole*Kse)) > 0) == 0
    disp('******** The closed-loop system is stable. ********')
end












