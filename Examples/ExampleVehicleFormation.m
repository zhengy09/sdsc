
%% *** Controller synthesis for networked system via chordal decompostion ******
% input: Gp --> graph for dynamical coupling; Gc --> communication graph
%        A{ij},B{i} --> dynamic matrices; 
% output: structured feedback gains K which corresponds to Gc;
%% ******************************************************************************
%%
clear; close all; clc

%% Data Generation Part
Node = 10:20:50;    %% the number of nodes in the graph
n = 2;     %% the dimension of state in each node
m = 1;     %% the dimension of input in each node

Bandwidth = 2;
Num = length(Node);
%% Time
TimeData = zeros(Num,2);
TimeSDP = zeros(Num,2);
TimeTotal = zeros(Num,2);
TimeGraph = zeros(Num,1);

CountFlag1 = 0;
CountFlag2 = 0;
for i = 1:length(Node)
    N = Node(i);
    Gp = zeros(N);     % decoupled dynamics
    Gc = zeros(N);     % banded communication
    for ii = 1:N
        for jj = 1:N
            if abs(ii -jj) <= Bandwidth && ii ~= jj
                Gc(ii,jj) = 1;Gc(jj,ii) = 1;
            end
        end
    end

    % Dynamic matrices
    A = cell(N,N);   % matrices for A
    B = cell(N,1);   % matrices for B

    % Dynamical part
    for ii = 1:N
        B{ii} = [0;1];
        A{ii,ii} = [0 1; 0 0];
    end

    %  load dataFailure
    [Awhole, Bwhole] = NetStateModel(A,B,Gp);
    % rank(ctrb(Awhole,Bwhole)) % this function is not accurate for a lagre N

    %% Design of structred feedback gains I: centralized solution (block diagonal Q)
    [K,Q,R,SedumiInfo,Tdata,Tsolver] = sgCentr(A,B,Gp,Gc);
    TimeData(i,1) = Tdata;
    TimeSDP(i,1) = Tsolver;
    TimeTotal(i,1) = Tdata+Tsolver;

    Flag = StrucCheck(K,Gc,m,n);
    if sum(real(eig(Awhole+Bwhole*K)) > 0) == 0
        disp('******** The closed-loop system is stable. ********')
        CountFlag1 = CountFlag1+1
    end


    %% Design of structred feedback gains III: Sequential design method
    %% Note that this will result undesesirable large gains
    [Kse,Qse,Rse,Tdatase,Tsolverse,Tgraphse,Ttotalse] = sgSeq(A,B,Gp,Gc);
    TimeData(i,2) = Tdatase;
    TimeSDP(i,2) = Tsolverse;
    TimeTotal(i,2) = Tdatase+Tsolverse;
    TimeGraph(i) = Tgraphse;

    Flag = StrucCheck(Kse,Gc,m,n);
    if sum(real(eig(Awhole+Bwhole*Kse)) > 0) == 0
        disp('******** The closed-loop system is stable. ********')
        CountFlag2 = CountFlag2 + 1
    end
end

