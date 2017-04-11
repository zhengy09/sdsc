function [LayNum, LayIndex] = LayerReform(Adj)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[n1,n1] = size(Adj);
% Atemp = triu(Adj);
Atemp = Adj;
% % LayNum = find(Atemp' ~= 0);
% LayIndex = cell(LayNum,1);

LayIndex{1} = 1;

%% calcualte the number already done
LayNum = length(LayIndex);
NodeNum = 0;
for i = 1:LayNum
    NodeNum = NodeNum + length(LayIndex{i});
end

while NodeNum < n1
    UpperLay = LayIndex{LayNum};  %% the node in upper layer
    Num = length(UpperLay);     %% the number of nodes in this layer
    for j = 1:Num
        Index1 = find(Atemp(UpperLay(j),:) == 1);
        Index = Index1;
        if LayNum > 1
            Index = setdiff(Index1,LayIndex{LayNum-1});
        end
        if length(LayIndex) == LayNum
            LayIndex{LayNum+1} = Index;
        else
            LayIndex{LayNum+1} = union(LayIndex{LayNum+1},Index);
        end
    end
    
    LayNum = length(LayIndex);
    NodeNum = 0;
    for i = 1:LayNum
        NodeNum = NodeNum + length(LayIndex{i});
    end
    if NodeNum == n1
        break;
    end

end

end