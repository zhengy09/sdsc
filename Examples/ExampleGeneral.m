
%% *** Controller synthesis for networked system via chordal decompostion ******
% input: Gp --> graph for dynamical coupling; Gc --> communication graph
%        A{ij},B{i} --> dynamic matrices; 
% output: structured feedback gains K which corresponds to Gc;
%% ******************************************************************************
%%
clear; close all; clc

%% Data Generation Part
n = 2;     %% the dimension of state in each node
m = 1;     %% the dimension of input in each node

% Chordal graph generation
Node = 10:50:200;    %% the number of nodes in the graph
Thresh = 5; 
Num = length(Node);

%% Time
TimeData = zeros(Num,2);
TimeSDP = zeros(Num,2);
TimeTotal = zeros(Num,2);
TimeGraph = zeros(Num,1);

CountFlag1 = zeros(Num,1);
CountFlag2 = zeros(Num,1);

p = 0.1; %% move 10% edges

for i = 1:Num
    
    N = Node(i);
    G = chordalGen(N,Thresh) - eye(N);       

    % Graph: Here we assume plant and communication share the same chordal graph
    Gp = GraphMoveEdge(G,p);     % graph for plant
    Gc = GraphAverageEdge(G,Gp);  % -- to make Gc cover Gp

   
    % Dynamic matrices
    A = cell(N,N);   % matrices for A
    B = cell(N,1);   % matrices for B

    % Dynamical part
    for ii = 1:N
        for jj = 1:N
            if Gp(ii,jj) == 1
                A{ii,jj} = exp(-(ii-jj)^2/10)*eye(n);     %% node j has influence on node i
            end  
            if ii == jj
                B{ii} = [0;1];
                A{ii,ii} = [1 1; 1 2];
            end      
        end
    end

    %  load dataFailure
    [Awhole, Bwhole] = NetStateModel(A,B,Gp);
   
    %% Design of structred feedback gains I: centralized solution (block diagonal Q)
    [K,Q,R,SedumiInfo,Tdata,Tsolver] = scCentr(A,B,Gp,Gc);
    TimeData(i,1) = Tdata;
    TimeSDP(i,1) = Tsolver;
    TimeTotal(i,1) = Tdata+Tsolver;

    Flag = StrucCheck(K,Gc,m,n);
    if sum(real(eig(Awhole+Bwhole*K)) > 0) == 0
        disp('******** The closed-loop system is stable. ********')
        CountFlag1(i) = 1;
    end

    %% Design of structred feedback gains III: Sequential design method
    % Note that this will result undesesirable large gains
    [Kse,Qse,Rse,Tdatase,Tsolverse,Tgraphse,Ttotalse] = scSeq(A,B,Gp,Gc);
    TimeData(i,2) = Tdatase;
    TimeSDP(i,2) = Tsolverse;
    TimeTotal(i,2) = Tdatase+Tsolverse;
    TimeGraph(i) = Tgraphse;

    Flag = StrucCheck(Kse,Gc,m,n);
    if sum(real(eig(Awhole+Bwhole*Kse)) > 0) == 0
        disp('******** The closed-loop system is stable. ********')
        CountFlag2(i) = 1;
    end
end
