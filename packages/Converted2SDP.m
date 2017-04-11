function [SDP,Tdata] = Converted2SDP(A,B,Gp,Gc)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here


[N,temp] = size(Gp);     % Number of nodes in the graph
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

% Ht = sparse(Nc,(N*n*2)^2);

Indexi = zeros((N*n*2)*10,1);        %% initial allocation
Indexj = zeros((N*n*2)*10,1);
Valueij = zeros((N*n*2)*10,1);

% contructing the basis of matrix               
[u,v] = find(triu(Qindex) == 1);
[u1,v1] = find(Rindex == 1);
Hi = zeros(N*n*2);
IndexNum = 0;
for i = 1:Nc
    if i < Nc1+1
        Q = zeros(N*n);                       % basis of Matrix Q
        Q(u(i),v(i)) = 1; Q(v(i),u(i)) = 1;
        
        tmp = Q*Awhole';
        Hi = blkdiag(-Q,tmp+tmp');
%         Hi(1:N*n,1:N*n) = -Q;      
%         Hi(N*n+1:end,N*n+1:end) = tmp+tmp';
        
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
%           Ht(i,:) = sparse(y'); 
        
    else
        R = zeros(N*m,N*n);                   % basis of Matrix R
        R(u1(i-Nc1),v1(i-Nc1)) = 1; 
        tmp = Bwhole*R;
        Hi = blkdiag(zeros(N*n),tmp+tmp');
        %Hi(1:N*n,1:N*n) = 0;
        
        %Hi(N*n+1:end,N*n+1:end) = tmp+tmp';
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
        
%         [jj] = find(y ~= 0); ii = i*ones(length(jj),1);
%         Indexi = [Indexi;ii];Indexj = [Indexj;jj'];
%         Valueij = [Valueij;y(jj)'];
%           Ht(i,:) = sparse(y');    
    end
end


IndexNz = find(Indexi~=0);
Ht = sparse(Indexi(IndexNz),Indexj(IndexNz),Valueij(IndexNz));   %% this operation has faster speed

%   sum(sum(full(Ht1-Ht)))
%% Controller design  -->   using Sedumi
epsilon = 0.001;
dSDP = N*n*2;      % the dimension of SDP
H0 = - epsilon*eye(dSDP); c = vec(H0); c = sparse(c);
b = zeros(Nc,1);          K.s = dSDP; pars.eps =1.0e-08;

SDP.At = Ht;
SDP.b = b;
SDP.c = c;
SDP.K = K;

Tdata = toc;


end

