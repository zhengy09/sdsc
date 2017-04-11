
%% Example 1: Stablization of dynamical systems over cyclic tree

clc;clear

n = 2;
N = 1:5;

TimeData = zeros(length(N),2);
TimeSDP = zeros(length(N),2);
TimeTotal = zeros(length(N),2);
TimeGraph = zeros(length(N),1);

for N = 2:3
    %% generate communication graph && banded communication
    Gc = zeros(sum(power(n,0:N)));      %% communication graph
    Dis = zeros(sum(power(n,0:N)),2);   %% positon of each node
    R = 20; Dis(1,1) = 0;Dis(1,2) = 0;
    Index = 1;
    for i = 1:N
        Angle = -pi/2+pi/(n^i):2*pi/(n^i):pi/2*3-pi/(n^i);
        for j = 1:n^i
            Dis(j+Index,1) = R*i*cos(Angle(j));Dis(j+Index,2) = R*i*sin(Angle(j));
        end   
        Index = Index + n^i;   
    end

    % Generate communication graph Gc
    Index = 0;
    for i = 0:N-1
        for j = 1:n^i
            CandiNode = Dis(Index+n^i+1:Index+n^i+n^(i+1),:);
            RDis = zeros(length(CandiNode),1);                     %% relative distance
            for ii = 1: length(CandiNode)
                RDis(ii) = norm(CandiNode(ii,:)- Dis(j+Index,:));
            end
            [Temp,PIndex] = sort(RDis);                             %% the neareast n node
            for ii = 1:n
                Gc(j+Index,PIndex(ii)+Index+n^i) = 1;Gc(PIndex(ii)+Index+n^i,j+Index) = 1;
            end
        end   
        Index = Index + n^i;   
    end
    
    Gc = triu(Gc,0);
    Gp = Gc';
      
    %% Dynamic matrices
    N1 = sum(power(n,0:N));   %% number of nodes;
    m = 1; n1 = 2;
    A = cell(N1,N1);   % matrices for A
    B = cell(N1,1);   % matrices for B

    % Dynamical part
    for i = 1:N1
        for j = 1:N1
            if Gp(i,j) == 1
                A{i,j} = exp(-(i-j)^2/10)*eye(n1);     %% node j has influence on node i
            end  
            if i == j
                B{i} = [0;1];
                A{i,i} = [1 1; 1 2];
            end      
        end
    end

    [Awhole, Bwhole] = NetStateModel(A,B,Gp);
    %% Design of structred feedback gains I: centralized solution (block diagonal Q)
    [K,Q,R,SedumiInfo,Tdata,Tsolver] = scCentr(A,B,Gp,Gc);
    Flag = StrucCheck(K,Gc,m,n1);
    if sum(real(eig(Awhole+Bwhole*K)) > 0) == 0
        disp('******** The closed-loop system is stable. ********')
    end

    TimeData(N,1) = Tdata;
    TimeSDP(N,1) = Tsolver;
    TimeTotal(N,1) = Tdata+Tsolver;

    %% Design of structred feedback gains II: Sequential design method
    % Note that this will result undesesirable large gains
    [Kse,Qse,Rse,Tdatase,Tsolverse,Tgraph,Ttotal] = scSeq(A,B,Gp,Gc);

    Flag = StrucCheck(Kse,Gc,m,n1);
    if sum(real(eig(Awhole+Bwhole*Kse)) > 0) == 0
        disp('******** The closed-loop system is stable. ********')
    end
    TimeData(N,2) = Tdatase;
    TimeSDP(N,2) = Tsolverse;
    TimeTotal(N,2) = Tdatase+Tsolverse;
    TimeGraph(N,1) = Tgraph;
end



