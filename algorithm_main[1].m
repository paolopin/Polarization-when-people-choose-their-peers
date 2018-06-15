% Matlab corde for
% ``Polarization when people choose their peers''
% Ugo Bolletta, Paolo Pin
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program allows to perform numerical simulation for the model, after
% having selected a value (or a range of values) for the parameters V,f,n, 
% and with a uniform distribution of initial opinions. The latter can be
% changed to any type of discrete distribution. 
% The algorithm selects the payoff maximizing network in time t, given a
% distribution of opinions at t-1. Thus, opinions are updated. The outcome
% is a vector of opinions with size nx1. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ADDITIONAL FILES NEEDED:
% - net_sel.m (the algorithm to create the optimal set of links for an 
%   agent)
% - networkComponents.m (used to count the number of components at the end 
%   of the simulation)
% - diameter.m (calculates the diameter of the network - SOURCE: Matlab 
%   Tools for Network Analysis (2006-2011)Copyright (c) 2011,
%   Massachusetts Institute of Technology.)
% - simple_dijkstra.m (used in diameter.m- SOURCE: Matlab Tools for Network
%   Analysis (2006-2011)Copyright (c) 2011, Massachusetts Institute of 
%   Technology.)
% - figure6.m to replicate the figure 6 in the paper
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTES: 
% In the ``Parameters'' section we propose several combination of
% parameters that allow to reproduce the figures in the paper. The
% initially uncommented section is left to the user to put in any
% combination of parameters.
clear 

%% Parameters %%
%In this section the reader can manipulate the main parameters to play
%around with the model. For replication purposes leave this section
%commented and uncomment further sections.

rng=1;
n=101; % number of nodes
T=50; % iterations
A=zeros(n,n);
f=0.5; % flexibility
V=0.07; % value of interactions
grid_fine=100;

%Figure 2a --- UNCOMMENT TO REPLICATE%
% rng=1;
% n=101; %number of nodes
% T=25; % iterations
% A=zeros(n,n);
% f=0.5;
% V=0.015;
% grid_fine=100;

%Figure 2b --- UNCOMMENT TO REPLICATE%
% rng=1;
% n=101; %number of nodes
% T=25; % iterations
% A=zeros(n,n);
% f=0.5;
% V=0.02;
% grid_fine=100;

%Figure 3a --- UNCOMMENT TO REPLICATE%
% rng=1;
% n=101; %number of nodes
% T=25; % iterations
% A=zeros(n,n);
% f=0.5;
% V=0.0350916;
% grid_fine=100;

%Figure 3b --- UNCOMMENT TO REPLICATE%
% rng=1;
% n=101; %number of nodes
% T=25; % iterations
% A=zeros(n,n);
% f=0.5;
% V=0.0350917;
% grid_fine=100;

%Figure 4a --- UNCOMMENT TO REPLICATE%
% rng=1;
% n=101; %number of nodes
% T=25; % iterations
% A=zeros(n,n);
% f=0.5;
% V=0.069;
% grid_fine=100;


%Figure 4b --- UNCOMMENT TO REPLICATE%
% rng=1;
% n=101; %number of nodes
% T=25; % iterations
% A=zeros(n,n);
% f=0.5;
% V=0.07;
% grid_fine=100;

%Figure 5 --- UNCOMMENT TO REPLICATE%
% n=80;
% f=linspace(0.01,0.95,grid_fine); % flexibility
% V=linspace(0.005,0.1,grid_fine); % benefit
% Pol_1=zeros(length(V)*length(f),T);
% Pol_2=zeros(length(V)*length(f),T);
 %V=(16*f-12*f^2-40*f*n+24*f^2*n+20*f*n^2-9*f^2*n^2)/(48*n^2)-0.01; %gross value of a link
%  V=(1/n)^2*(f*(1-f))/(1+f)^2;
S=zeros(n,T); %size of peer group
M=zeros(T,n);
C=zeros(n,T);
y=zeros(n,n);
B=zeros(n,T);
Z=zeros(length(f),length(V));
W=zeros(length(f),length(V));
uW=zeros(length(f),length(V));
Y=zeros(length(f),length(V));
grid=zeros(length(f), length(V)); 
for i=1:n
    M(1,i)=(i-1)/(n-1); % Creates a uniform distribution of initial opinions
end

%% Network formation

for w=1:length(f)
for z=1:length(V)
for t=2:T
    a=M(t-1,:); % (vector of opinions)
    for i=1:n
       
        F=net_sel(a,i,n,f(w),V(z));
        h=sum(F>0);
        S(i,t)=h;
        if h>0
            m=mean(a(logical(F)));
        end
        M(t,i)=(1-f(w))*a(i)+f(w)*m;
        A(i,:)=F;
     end
                if t==2
                    [Z(w,z),sizes,member] = networkComponents(A);
                     W(w,z)=diameter(A);
                     uW(w,z)=diameter(((A'+A)>0)+0);
                end
                if t==T
                    Y(w,z)=diameter(A);
                end
 B(:,t)=M(t,:)-M(t-1,:);
end
 x=sign(B);
 for i=1:n
     for j=1:T-1
        y(i,j)=x(i,j)*x(i,j+1);
     end
 end
 grid(w,z)=sum(sum((y(:)<0)));
 if  any(y(:)<0)==1
    grid(w,z)=1;
 else
    grid(w,z)=0;
 end
    figure
    for ii=1:n
    plot(M(:,ii));
    hold on
    end
end 

datetime('now')
w

end

save('results')



