function [K, Q, R ] = scCentrYalmip(A,B,Gp,Gc,alpha,beta,gamma)
% Design of Structred Feedback Gains via centralized solution
% Input: two graphs: Gp, graph for plant; Gc, graph for communication 
%        Matrices for Dynamic, A,B (cell form)
% Assumption: Block Digonal Q


[N,temp] = size(Gp);     % Number of nodes in the graph
[n,m] = size(B{1});      % dimensions of input and state in each node 

%% Reformulation of data --> the whole state space model & chordal graph Gex
[Awhole, Bwhole] = NetStateModel(A,B,Gp);

tic
%% Controller design  --> sparsity pattern in Q,R

Q = sdpvar(n);               %% block diagonal Q 
for i = 2:N
    Q = blkdiag(Q,sdpvar(n));
end
R = sdpvar(N*m,N*n);        %% Matrix R has sparsity pattern in Gc
for i = 1:N                     
    for j = 1:N
        if Gc(i,j) ~= 1 && i ~= j
            R((i-1)*m+1:i*m,(j-1)*n+1:j*n) = zeros(m,n);
        else
            %R((i-1)*m+1:i*m,[(j-1)*n+2,j*n]) = 0;  %% part of state feedback
        end
    end
end

epsilon = 0.001;

% alpha = 0.01;
% 
% beta  = 100;
% gamma = 10;

opts = sdpsettings('verbose',1,'solver','sedumi');
%% centralized solution
Const = [Q - 1/gamma*eye(n*N)>= 0];
Const = [Const, Q*Awhole'+Awhole*Q+R'*Bwhole'+Bwhole*R + epsilon*eye(n*N) + alpha*Q<=0];

for i = 1:N                     
    for j = 1:N
        if Gc(i,j) == 1 || i == j
           Const =[Const, [-beta*eye(n) R((i-1)*m+1:i*m,(j-1)*n+1:j*n)'; R((i-1)*m+1:i*m,(j-1)*n+1:j*n) -eye(m)] <=0];
        end
    end
end


Obj = 0;
sol = optimize(Const,Obj,opts);

Q = value(Q);
R = value(R);
K = R*Q^(-1);

end

