function [K, Q, R ] = scSeqYalmip(A,B,Gp,Gc,alpha,beta,gamma)
% Design of Structred Feedback Gains via sequential design
% Input: two graphs: Gp, graph for plant; Gc, graph for communication 
%        Matrices for Dynamic, A,B (cell form)
% Assumption: Block Digonal Q


[N,temp] = size(Gp);     % Number of nodes in the graph
[n,m] = size(B{1});      % dimensions of input and state in each node 

%% Reformulation of data --> the whole state space model & chordal graph Gex
[Awhole, Bwhole] = NetStateModel(A,B,Gp);

Q = sdpvar(n);               %% block diagonal Q 
for i = 2:N
    Q = blkdiag(Q,sdpvar(n));
end
R = sdpvar(N*m,N*n);        %% Matrix R has sparsity pattern in Gc
for i = 1:N                     
    for j = 1:N
        if Gc(i,j) ~= 1 && i ~= j
            R((i-1)*m+1:i*m,(j-1)*n+1:j*n) = zeros(m,n);
        end
    end
end

epsilon = 0.001;
%% Clique 1

Q1 = Q(1:2*n,1:2*n);
R1 = R(1:2*m,1:2*n); 
Factor1 = ones(n*2);
Factor1(n+1:2*n,n+1:2*n) = 1/2;

opts = sdpsettings('verbose',1,'solver','sedumi');
Const = [Q1 - 1/gamma*eye(n*2)>= 0];
Const = [Const, (Q1*Awhole(1:2*n,1:2*n)'+Awhole(1:2*n,1:2*n)*Q1 + ...
           R1'*Bwhole(1:2*n,1:2*m)'+Bwhole(1:2*n,1:2*m)*R1 + epsilon*eye(n*2) + alpha*Q1).*Factor1<=0];

for i = 1:2                     
    for j = 1:2
        if Gc(i,j) == 1 || i == j
           Const =[Const, [-beta*eye(n) R1((i-1)*m+1:i*m,(j-1)*n+1:j*n)'; R1((i-1)*m+1:i*m,(j-1)*n+1:j*n) -eye(m)] <=0];
        end
    end
end

Obj = 0;
sol = optimize(Const,Obj,opts);

Q(1:2*n,1:2*n) = value(Q1);
R(1:2*m,1:2*n) = value(R1);

%% Clique 2
Q2 = Q(n+1:3*n,n+1:3*n);
R2 = R(m+1:3*m,n+1:3*n); 
Factor2 = ones(n*2);
Factor2(1:n,1:n) = 1/2;

opts = sdpsettings('verbose',1,'solver','sedumi');
Const = [Q2 - 1/gamma*eye(n*2)>= 0];
Const = [Const, (Q2*Awhole(n+1:3*n,n+1:3*n)'+Awhole(n+1:3*n,n+1:3*n)*Q2 + ...
           R2'*Bwhole(n+1:3*n,m+1:3*m)'+Bwhole(n+1:3*n,m+1:3*m)*R2 + epsilon*eye(n*2) + alpha*Q2).*Factor2<=0];

for i = 2:3                     
    for j = 2:3
        if (Gc(i,j) == 1 || i == j)
            if i == 2 && j == 2 
                
            else 
                Const =[Const, [-beta*eye(n) R2((i-2)*m+1:(i-1)*m,(j-2)*n+1:(j-1)*n)'; R2((i-2)*m+1:(i-1)*m,(j-2)*n+1:(j-1)*n) -eye(m)] <=0];
            end
         end
    end
end

Obj = 0;
sol = optimize(Const,Obj,opts);

Q(n+1:3*n,n+1:3*n) = value(Q2);
R(m+1:3*m,n+1:3*n) = value(R2);

Q = value(Q);
R = value(R);
K = R*Q^(-1);

end

