function [K, Q, R, SedumiInfo,Tdata, Tsolver] = scCentr(A,B,Gp,Gc,Opts)
% Design of Structred Feedback Gains via centralized solution
% Input: two graphs: Gp, graph for plant; Gc, graph for communication 
%        Matrices for Dynamic, A,B (cell form)
%        Opts: 1 using SeDuMi; 2 using SparseCoLo
% Assumption: Block Digonal Q

if(nargin <= 4)
    Opts = 1;
end

[N,~] = size(Gp);     % Number of nodes in the graph
[n,m] = size(B{1});      % dimensions of input and state in each node 

%% Reformulation of data --> the whole state space model & chordal graph Gex
[Awhole, Bwhole] = NetStateModel(A,B,Gp);

tic
%% Controller design  --> sparsity pattern in Q,R
Qindex = ones(n);               %% block diagonal Q 
for i = 2:N
    Qindex = blkdiag(Qindex,ones(n));
end
Rindex = zeros(N*m,N*n);        %% Matrix R has sparsity pattern in Gc
for i = 1:N                     
    for j = 1:N
        if Gc(i,j) == 1 || i == j
            Rindex((i-1)*m+1:i*m,(j-1)*n+1:j*n) = ones(m,n);
        end
    end
end

%% Controller design  --> reformulate the data into standard SDP form in Sedumi
Nc1 = sum(sum(triu(Qindex)));       % The number of constraints for Q
Nc = Nc1 + sum(sum(Rindex));        % The number of constraints for R

Indexi = zeros((N*n*2)*10,1);        %% initial allocation
Indexj = zeros((N*n*2)*10,1);
Valueij = zeros((N*n*2)*10,1);

% contructing the basis of matrix               
[u,v] = find(triu(Qindex) == 1);
[u1,v1] = find(Rindex == 1);

IndexNum = 0;
for i = 1:Nc
    if i < Nc1+1
        Q = zeros(N*n);                       % basis of Matrix Q
        Q(u(i),v(i)) = 1; Q(v(i),u(i)) = 1;
        Hi = blkdiag(-Q,Q*Awhole'+Awhole*Q);
        y = vec(Hi);
              
        [jj] = find(y ~= 0); ii = i*ones(length(jj),1);
        
        if IndexNum+length(jj)+1 > length(Indexi)
            Indexi = [Indexi;zeros((N*n*2)*10,1)];
            Indexj = [Indexj;zeros((N*n*2)*10,1)];
            Valueij = [Valueij;zeros((N*n*2)*10,1)];
        end
        
        Indexi(IndexNum+1:IndexNum+length(jj)) = ii;
        Indexj(IndexNum+1:IndexNum+length(jj)) = jj;
        Valueij(IndexNum+1:IndexNum+length(jj)) = y(jj);
        IndexNum = IndexNum + length(jj);        
    else
        R = zeros(N*m,N*n);                   % basis of Matrix R
        R(u1(i-Nc1),v1(i-Nc1)) = 1; 
        Hi = blkdiag(zeros(N*n),Bwhole*R+R'*Bwhole');
        y = vec(Hi)';
        
        [jj] = find(y ~= 0); ii = i*ones(length(jj),1);
        
        if IndexNum+length(jj)+1 > length(Indexi)
            Indexi = [Indexi;zeros((N*n*2)*10,1)];
            Indexj = [Indexj;zeros((N*n*2)*10,1)];
            Valueij = [Valueij;zeros((N*n*2)*10,1)];
        end
        
        Indexi(IndexNum+1:IndexNum+length(jj)) = ii;
        Indexj(IndexNum+1:IndexNum+length(jj)) = jj;
        Valueij(IndexNum+1:IndexNum+length(jj)) = y(jj);
        IndexNum = IndexNum + length(jj); 
    end
end


IndexNz = find(Indexi~=0);
Ht = sparse(Indexi(IndexNz),Indexj(IndexNz),Valueij(IndexNz));   %% this operation has faster speed

%   sum(sum(full(Ht1-Ht)))
%% Controller design  -->   using Sedumi
epsilon = 0.001;
dSDP = N*n*2;      % the dimension of SDP
H0 = - epsilon*eye(dSDP); c = vec(H0); c = sparse(c);
b = zeros(Nc,1);          K.s = dSDP; 
Tdata = toc;

if Opts == 1  %% solving by SeDuMi
    tic
    pars.eps =1.0e-08;
    [x,y,SedumiInfo] = sedumi(Ht,b,c,K,pars);
    Tsolver = toc;
else      %% solving by SparseCoLo
    tic
    K.s = dSDP;
    J.f = Nc; 
    [x,y,infoCoLO,cliqueDomain,cliqueRange,LOP] = sparseCoLO(Ht,b,c,K,J);
    Tsolver = toc;
end
%% Controller design  -->   converting into controller form
[Q,R] = Vector2Matrix(y,u,v,u1,v1,N,n,m,Nc,Nc1);     
K = R*Q^(-1);

end

