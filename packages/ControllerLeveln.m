function [Rs,Qs,Tdata,Tsolver] = ControllerLeveln(MC1,MCindex,A,B,Gp,Gc,Rs,Qs,NodeTime,epsilon)
%   Subfunction for sovling the controller gains in each clique
%   Input: maximal cliques, MC1; Dynamic matrices: A,B;
%          graphs,GP,Gc;    Controller gains for computing, Qd,Rd;
%          Overlapping elements, NodeTime

    [n,m] = size(B{1});
    NodeIndex = find(MC1(:,MCindex) == 1);  %% The Index of nodes in Level 0, Rij, Qij, i,j --> NodeIndex
    NodeNum = length(NodeIndex);            %% the number of nodes in this clique

    Qindex = zeros(n*NodeNum); Rindex = zeros(m*NodeNum,n*NodeNum);
    Qdone = zeros(n*NodeNum); Rdone = zeros(m*NodeNum,n*NodeNum);
    Atemp = zeros(n*NodeNum); Btemp = zeros(n*NodeNum,m*NodeNum);

    %% Dynamcis & Controller structure
    for ii = 1:length(NodeIndex)   %% Dynamics
        Btemp((ii-1)*n+1:ii*n,(ii-1)*m+1:ii*m) = B{NodeIndex(ii)}; 
        for jj = 1:length(NodeIndex)
            if Gp(NodeIndex(ii),NodeIndex(jj)) == 1 || ii == jj
                Atemp((ii-1)*n+1:ii*n,(jj-1)*n+1:jj*n) = A{NodeIndex(ii),NodeIndex(jj)};
            end
        end
    end
    
    % controller structure
    for ii = 1:length(NodeIndex)   %% Dynamics
        if isempty(Qs{NodeIndex(ii)})    %% this variable is free
            Qindex((ii-1)*n+1:ii*n,(ii-1)*n+1:ii*n) = ones(n);
        else
            Qdone((ii-1)*n+1:ii*n,(ii-1)*n+1:ii*n) = Qs{NodeIndex(ii)};
        end
        for jj = 1:length(NodeIndex)
            if Gc(NodeIndex(ii),NodeIndex(jj)) == 1 || ii == jj
                if isempty(Rs{NodeIndex(ii),NodeIndex(jj)}) 
                    Rindex((ii-1)*m+1:ii*m,(jj-1)*n+1:jj*n) = ones(m,n);
                else
                    Rdone((ii-1)*m+1:ii*m,(jj-1)*n+1:jj*n) = Rs{NodeIndex(ii),NodeIndex(jj)};
                end
            end
        end
    end
 
 %% converting into standard form
    tic
    Nc1 = sum(sum(triu(Qindex)));       % The number of constraints for Q
    Nc = Nc1 + sum(sum(Rindex));        % The number of constraints for R
    Ht = sparse(Nc,(NodeNum*n*2)^2);
    
    % Equally divide the overlapping elements
    Scale =blkdiag(ones(n*NodeNum),kron(1./NodeTime(NodeIndex,NodeIndex),ones(n)));          
    
    [u,v] = find(triu(Qindex)); [u1,v1] = find(Rindex);
    for ii = 1: Nc
        if ii < Nc1+1
            Qtemp = zeros(n*NodeNum);
            Qtemp(u(ii),v(ii)) = 1; Qtemp(v(ii),u(ii)) = 1;
            tmp = Qtemp*Atemp';
            Hi = blkdiag(-Qtemp,tmp+tmp');
            Hi = Scale.*Hi;
            y = vec(Hi)';
            Ht(ii,:) = sparse(y);
        else
            Rtemp = zeros(m*NodeNum,n*NodeNum);
            Rtemp(u1(ii-Nc1),v1(ii-Nc1)) = 1;
            tmp=Rtemp'*Btemp';
            Hi = blkdiag(zeros(NodeNum*n),tmp+tmp');
            Hi = Scale.*Hi;
            y = vec(Hi)';
            Ht(ii,:) = sparse(y);
        end
    end


    %% using Sedumi
    dSDP = NodeNum*n*2;      % the dimension of SDP
    Temp1 = blkdiag(-Qdone,Qdone*Atemp'+Atemp*Qdone);
    Temp2 = blkdiag(zeros(NodeNum*n),Rdone'*Btemp'+Btemp*Rdone);
    H0 = - Scale.*(epsilon*eye(dSDP) + Temp1 + Temp2);
    c = vec(H0); c = sparse(c);
    b = zeros(m,1);
    Tdata = toc;

    K.s = dSDP; pars.eps =1.0e-08;pars.fid = 1;
    tic
    [x,y,info] = sedumi(Ht,b,c,K,pars);
    Tsolver = toc;
    
%     if info.feasratio < -0.5
%         error('infeasible')
%     end

    %% %% Recovery for the Rij and Qij
    [Qtemp,Rtemp] = Vector2Matrix(y,u,v,u1,v1,NodeNum,n,m,Nc,Nc1);        
    for ii = 1:length(NodeIndex)   
        if isempty(Qs{NodeIndex(ii)})    %% this variable is free
            Qs{NodeIndex(ii)} = Qtemp((ii-1)*n+1:ii*n,(ii-1)*n+1:ii*n);
        end      
        for jj = 1:length(NodeIndex)
            if Gc(NodeIndex(ii),NodeIndex(jj)) == 1 || ii == jj
                if isempty(Rs{NodeIndex(ii),NodeIndex(jj)}) 
                    Rs{NodeIndex(ii),NodeIndex(jj)} = Rtemp((ii-1)*1+1:ii*1,(jj-1)*n+1:jj*n);
                end
            end
        end
    end  

end

