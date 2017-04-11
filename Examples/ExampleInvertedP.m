
%% *** Controller synthesis for networked system via chordal decompostion ******
% input: Gp --> graph for dynamical coupling; Gc --> communication graph
%        A{ij},B{i} --> dynamic matrices; 
% output: structured feedback gains K which corresponds to Gc;
%% ******************************************************************************
%%
clear; close all;% clc

DataInverted

% communication graph
% case 1: decentralized
Gc = [1 1 0;
      1 1 1;
      0 1 1];

[Awhole, Bwhole] = NetStateModel(A,B,Gp);

%% Time
% TimeData = zeros(1,2); TimeSDP = zeros(1,2); TimeTotal = zeros(1,2); TimeGraph = zeros(1,1);

alpha = 0.5;
beta  = 10;
gamma = 10;

%% Design of structred feedback gains I: centralized solution (block diagonal Q)
[K, Q, R ] = scCentrYalmip(A,B,Gp,Gc,alpha,beta,gamma);

eig(Awhole+Bwhole*K)


%% Design of structred feedback gains III: Sequential design method
%% Note that this will result undesesirable large gains
[Kse, Qse, Rse ] = scSeqYalmip(A,B,Gp,Gc,alpha,beta,gamma);

eig(Awhole+Bwhole*Kse)


%%
t = 0:0.1:19.5; close all

%alpha = 1;

% x0 = [0.05; 0; 0.2; 0;0.05; 0; 0.2; 0;0.05; 0; 0.2; 0;];
x0 = [0.05; 0; 0.2; 0;0.0; 0; 0; 0;0; 0; 0; 0;];
x1 = zeros(length(t),length(x0));
x1n = zeros(length(t),1);
x2 = zeros(size(x1));
x2n = zeros(length(t),1);
x3 = zeros(length(t),1);
for i = 1:length(t)
    x1(i,:) = expm((Awhole+Bwhole*Kse)*t(i))*x0;
    x2(i,:) = expm((Awhole+Bwhole*K)*t(i))*x0;
    x1n(i) = norm(x1(i,:),'fro'); 
    x2n(i) = norm(x2(i,:),'fro');
    x3(i) = 2.5*exp(-alpha/2*t(i))* norm(x0,'fro');
end


figure; h1 = plot(t,x2n,'m','linewidth',1); hold on
h2 = plot(t,x1n,'b','linewidth',1); 
h3 = plot(t,x3,'k','linewidth',1); 

h = legend([h1,h2,h3],'Centralized','Sequentail','Guaranteed performance','Location','Northwest');
set(h,'box','off')
set(gcf,'Position',[250 150 350 350]);box off, grid off; %xlim([10 200]); ylim([1,10^4])
%set(gca,'XTick',[10,20,50,100,150,200]) 
xlabel('Time (s)');ylabel('Norm of error |x(t)|')
