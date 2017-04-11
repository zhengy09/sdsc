function [K,Q,R,Tdata,Tsolver,Tgraph,Ttotal] = scSeq(A,B,Gp,Gc)
%   Sequential design method based on the clique tree
%   Detailed explanation goes here

[N,temp] = size(Gp);     % Number of nodes in the graph
[n,m] = size(B{1});      % dimensions of input and state in each node 
Qs = cell(N,1); Rs = cell(N);   % the parameters in final gains

Tdata  = 0; Tsolver = 0;    %% for counting the consuming time
Tgraph = 0; Ttotal = [];   %% time for graph operation
%% RefoRrmulation of data --> chordal extension --> Gex
tic
Gs = Gp | Gp' | Gc | Gc';                        % super-graph: modelling connections in the whole system
Gex = chordalExt(Gs+eye(size(Gs)));              % chordal extension

[MC0,Adj,MC1] = prim(Gex);         % get the maximal cliques and cliue tree in the extended chordal graph
NodeTime = OverlappingTimes(Gex);  % The times of overlapping elements  -- operations in graph level

%% First to count the number of layers in the clique tree
% [Adj1,MC1] = TreeRegulation(Adj,MC1);    % re-index the maximal cliques
[LayNum,LayIndex] = LayerReform(Adj);          % get the index in each layer
Adj1 = Adj;            % by 20160104

%MC = maximalCliques(G);

% ELevelNum = sum(triu(Adj1)');            % The number of cliques in each level  --> Not exact, but this works
% LevelNum = find(ELevelNum ~= 0);         % The number of levels -1 

Tgraph = toc;
Tdata = Tgraph;
%%
epsilon = 0.01;
for i = 1:LayNum    % Breadth first,--> we first try to slove the controller for each level
    if i == 1                 %% Level 0 --> root node
        [Rs,Qs,Tdata0,Tsolver0] = ControllerLevel0(MC1,A,B,Gp,Gc,Qs,Rs,NodeTime,epsilon);
        Tdata = Tdata + Tdata0;
        Tsolver = Tsolver + Tsolver0;
        Ttotal = [Ttotal,Tsolver0];
    else                      %% other levels
        for ii = 1:length(LayIndex{i})   %% compute the controller for all cliques in this level
%             AdjTemp = triu(Adj1);
%             Index = find(AdjTemp(LevelNum(i),:) == 1);
            Index = LayIndex{i};
            MCindex = Index(ii);            %% Find the index of the clique
            %% solving this clique
            [Rs,Qs,Tdatan,Tsolvern] = ControllerLeveln(MC1,MCindex,A,B,Gp,Gc,Rs,Qs,NodeTime,epsilon);
            Tdata = Tdata + Tdata0;
            Tsolver = Tsolver + Tsolvern;
            Ttotal = [Ttotal,Tsolvern];
        end       
    end      
end

%% recover the data
R = zeros(m*N,N*n);Q = zeros(n*N);
for i = 1:N
    for j = 1:N
        if Gc(i,j) == 1 || i == j
            R((i-1)*1+1:i*1,(j-1)*n+1:j*n) = Rs{i,j};
        end
        if Gc(i,j) == 0 && i ~= j && ~isempty(Rs{i,j})
%             error('error')
        end
    end
    Q((i-1)*n+1:i*n,(i-1)*n+1:i*n) = Qs{i};
end

K = R*Q^(-1);

end

