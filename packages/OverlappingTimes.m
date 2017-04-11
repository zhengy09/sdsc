function [NodeTime] = OverlappingTimes( G )
%   UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

MC = maximalCliques(G); 
NodeTime = zeros(size(G));
[Temp,NumMc] = size(MC);
for i = 1:NumMc   %% calculate the times of each elements
    NodeIndex = find(MC(:,i) == 1);
    NodeNum = length(NodeIndex);
    for ii = 1:NodeNum
        for jj = ii:NodeNum
            if ii == jj
                NodeTime(NodeIndex(ii),NodeIndex(jj)) = NodeTime(NodeIndex(ii),NodeIndex(jj))+1;
            else
                NodeTime(NodeIndex(ii),NodeIndex(jj)) = NodeTime(NodeIndex(ii),NodeIndex(jj))+1;
                NodeTime(NodeIndex(jj),NodeIndex(ii)) = NodeTime(NodeIndex(jj),NodeIndex(ii))+1;
            end            
        end      
    end    
end

end

