function Flag = StrucCheck(K,Gc,m,n )
%   Check whether the controller is compatible with communication graph
%   Flag = 0, meanning the check is success

    [Node, temp] = size(Gc);
    Flag = 0;
    for i = 1:Node
        for j = 1: Node
            if Gc(i,j) == 0 & i ~= j & sum(sum(K((i-1)*m+1:i*m,(j-1)*n+1:j*n))) >0
                Flag = Flag + 1;
                warning('Something is wrong with Controller')
            end
        end
    end
    if Flag == 0
        disp(' ')
        disp('******** The structure of Controller K is compatible with communication graph Gc *********')
        disp(' ')
    end
end

