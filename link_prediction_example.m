close all
clear all
clc

addpath functions/
addpath tensor_toolbox/

% load("datasets/SmallW_main.mat","W"); dataset = 'SmallW';
load('datasets/UKfaculty.mat', 'W'); dataset='UKfaculty';

A = W;
G = graph(A);
G = max_connected_subgraph(G);
m = numedges(G);
n = numnodes(G);

figure, plot(G,'NodeLabel',{});


%% Set parameters 
c = .45;                    % pagerank teleport coeff
tau = .1;                   % remove tau% of edges
sigma = 1;                  % predict sigma% of removed edges


%% Number of random trials 
numtrials = 10;

%% Run test for a fixed value of alpha and p varying within parray
parray = [0 .25 .5 1 1.5 2];          
alpha0 = .5;
for j = 1 : numtrials
    ind_deleted_edges = randi([1,m],floor(tau*m),1);
    
    H = G.rmedge(ind_deleted_edges);
    A = H.adjacency();
    T = build_triangles_tensor(A,'type','random_walk');
    D = 1./sum(A,2);
    D(D == inf) = 0;
    D = spdiags(D,0,n,n);
    M = D*A;
    
    [score(j),~,~] = linear_pr_linkpredict(G,ind_deleted_edges,c,sigma);
    
    for jj = 1 : length(parray)
        p = parray(jj); 
        [score_nonlinear(jj,j),~,~] = nonlinear_pr_linkpredict(G,T,M',ind_deleted_edges,c,sigma,alpha0,p); 
    end    
    if (mod(j,5)==0 || j==1), fprintf('alpha fixed - trial number %d\n is over', j); end  
end

figure
boxplot((score_nonlinear ./ score)', 'Labels', parray,'Widths',.2); 
xlab = sprintf("alpha  = %1.1f, p", alpha0); xlabel(xlab);
ylabel('Relative performance ratio vs standard seeded PageRank');
yline(1,'--');


%% Run test for a fixed value of p and alpha varying within aarray
aarray = [.1 .2 .5 .8 .9 1];
p0 = 0;
for j = 1 : numtrials
    ind_deleted_edges = randi([1,m],floor(tau*m),1);
    
    H = G.rmedge(ind_deleted_edges);
    A = H.adjacency();
    T = build_triangles_tensor(A,'type', 'random_walk');
    D = 1./sum(A,2);
    D(D == inf) = 0;
    D = spdiags(D,0,n,n);
    M = D*A;
    
    for jj = 1 : length(aarray)
        alpha = aarray(jj); 
        [score_nonlinear_2(jj,j),~,~] = nonlinear_pr_linkpredict(G,T,M',ind_deleted_edges,c,sigma,alpha,p0);
    end    
    if (mod(j,5)==0 || j==1), fprintf('p fixed - trial number %d\n is over', j); end  
end

figure
boxplot((score_nonlinear_2 ./ score)', 'Labels', aarray,'widths',.2); 
xlab = sprintf("p  = %1.1f, alpha", p0); xlabel(xlab);
ylabel('Relative performance ratio vs standard seeded PageRank');
yline(1,'--');







