function [Flag] = ExampleGraphPlot(G,Dis,Msize)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%% Msize = 26

figure;
plot(Dis(:,1),Dis(:,2),'k.','markersize',Msize,'linewidth',1.5); hold on;
set(gcf,'position',[300,200,400,400])
set(gca,'position',[0.02,0.02,0.96,0.96])
axis off
box off;
xlim([min(Dis(:,1))-1,max(Dis(:,1))+1]);
ylim([min(Dis(:,2))-1,max(Dis(:,2))+1]);


%% plot links
[Num,temp] = size(G);
for i = 1: Num
    for j = i+1:Num
        if G(i,j) == 1
            plot([Dis(i,1),Dis(j,1)],[Dis(i,2),Dis(j,2)],'k','linewidth',1);hold on;
        end
    end    
end

% figure
% Kc = kron(G+eye(size(G)),ones(1,2));
% 
% [u,v] = find(Kc ~=0);
% plot(u,v,'m.','markersize',20,'linewidth',1.5);
% box off?
% set(gcf,'position',[300,200,400,600])


end

