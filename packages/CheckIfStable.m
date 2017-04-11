function [Flag,Awhole,Bwhole,Qmc,Rmc] = CheckIfStable(A,B,Gp,Gc,Nodedone,Q,R)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    NumNode = length(Nodedone);
    Amc = cell(NumNode);  Bmc = cell(NumNode,1);
    Gpmc = Gp(Nodedone,Nodedone);
    Gcmc = Gc(Nodedone,Nodedone);
    [n,m] = size(B{1});
    Qmc = zeros(NumNode*n);Rmc = zeros(NumNode*m,NumNode*n);
    for ii = 1:NumNode
        for jj = 1:NumNode
            if ii == jj
                Amc{ii,ii} = A{Nodedone(ii),Nodedone(ii)};
                Bmc{ii} = B{Nodedone(ii)};
            elseif Gpmc(ii,jj) == 1
                Amc{ii,jj} = A{Nodedone(ii),Nodedone(jj)};
            end
            Qmc((ii-1)*n+1:ii*n,(jj-1)*n+1:jj*n) = Q((Nodedone(ii)-1)*n+1:Nodedone(ii)*n,(Nodedone(jj)-1)*n+1:Nodedone(jj)*n);
            Rmc((ii-1)*m+1:ii*m,(jj-1)*n+1:jj*n) = R((Nodedone(ii)-1)*m+1:Nodedone(ii)*m,(Nodedone(jj)-1)*n+1:Nodedone(jj)*n);
        end
    end
    
    [Awhole, Bwhole] = NetStateModel(Amc,Bmc,Gpmc);
    if sum(eig(Awhole*Qmc+Qmc*Awhole'+Bwhole*Rmc+Rmc'*Bwhole')>0)
        warning('controller is not right!')
        Flag = 1;
    else
        Flag = 0;
    end
end

