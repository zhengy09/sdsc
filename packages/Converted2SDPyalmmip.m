function [SDP,Tdata] = Converted2SDPyalmmip(A,B,Gp,Gc)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here


[N,temp] = size(Gp);     % Number of nodes in the graph
[n,m] = size(B{1});      % dimensions of input and state in each node 

[Awhole, Bwhole] = NetStateModel(A,B,Gp);

tic
%--> sparsity pattern in Q,R
Q = sdpvar(n,n);               %% block diagonal Q 
for i = 2:N
    Q = blkdiag(Q,sdpvar(n,n));
end
R = sdpvar(m,n);        %% Matrix R has sparsity pattern in Gc
%Q = sdpvar(n,n);               %% block diagonal Q 
for i = 2:N
    R = blkdiag(R,sdpvar(m,n));
end
%R(1:end,1:end) = zeros(N*m,N*n);
for i = 1:N                     
    for j = 1:N
        if Gc(i,j) == 1
            R((i-1)*m+1:i*m,(j-1)*n+1:j*n) = sdpvar(m,n);
        else
           % R((i-1)*m+1:i*m,(j-1)*n+1:j*n) = zeros(m,n);
        end
    end
end

%% Controller design  --> reformulate the data into standard SDP form in Sedumi
epsilon = 0.001;
dSDP = N*n;      % the dimension of SDP
cnstr = [Q*Awhole'+Q*Awhole+R'*Bwhole' + Bwhole*R <=0,Q - epsilon*eye(dSDP)>=0];
obj = 0;
opts = sdpsettings('solver','sedumi');

[mod,~] = export(cnstr,obj,opts);

SDP.At = mod.A;
SDP.b = mod.b;
SDP.c = mod.C;
SDP.K = mod.K;

Tdata = toc;


end

