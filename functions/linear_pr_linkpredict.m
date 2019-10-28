function [score,added,X] = linear_pr_linkpredict(G,ind_deleted_edges,c,sigma)
%INPUT:
% G = the graph (digraph), 
% ind_deleted_edges = indices of edges to be removed edges, 
% c = the teleport constant for PR
% sigma = percentage of edges to predict


m = numedges(G);
n = numnodes(G);
e = ones(n,1);


H = G.rmedge(ind_deleted_edges);
% Standard pagerank
A = H.adjacency();
D = 1./(A*e);
D(D == inf) = 0;
D = spdiags(D,0,n,n);
P = D*A;


I = eye(n,n);
% Add in P rows of (1/n,...,1/n) on the empty rows so that P1 = 1
P(diag(D) == 0,:) = ones(sum(diag(D)==0),n)./n;
% Compute the similarity matrix by doing (1-c) (I - c*P)^{-1}
X = (1-c)*( (I-c*P')\I );
% We do not want to add loops in the graph and we dont want to consider
% edges that are already there
X(I>0) = -Inf;
X(A>0) = -Inf;
X = X + X';

% Take the k biggest elements in X(i,j) 
k = floor(sigma*length(ind_deleted_edges));
[~,J] = maxk(X(:),k);
[i,j] = ind2sub([n n],J);
added = [i,j];

score = ismember(added,G.Edges.EndNodes(ind_deleted_edges,:),'rows');
score = sum(score);

end

