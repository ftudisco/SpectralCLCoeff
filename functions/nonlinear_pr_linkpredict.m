function [score,added,X] = nonlinear_pr_linkpredict(G,T,M,ind_deleted_edges,c,sigma,alpha,p)
%INPUT:
% G = the graph (digraph), 
% ind_deleted_edges = indices of edges to be removed edges, 
% c = the teleport constant for PR
% sigma = percentage of edges to predict


% If alpha is equal to 1, the nonlinear diffusion map boils down to the
% standard linear pagerank matrix
if alpha == 1 
    [score,added,X] = linear_pr_linkpredict(G,ind_deleted_edges,c,sigma);
    return;
end


n = numnodes(G);
I = speye(n,n);
%% Compute the similarity vectors x_{k+1} = cM(x_k) + (1-c)v
maxiter = 10; 
X = zeros(n, n);
for i = 1 : n
   v = I(:,i);
   x = zeros(n,1);
   for k = 1 : maxiter
       x = c * ( alpha*M*x + (1-alpha)*Tp(T,x,p) ) + (1-c) * v;
       x = x ./ norm(x,1);
   end
   X(:,i) = x;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

X(I>0) = -Inf;
X(M>0) = -Inf;
X = X + X';

%% compute edges to add
k = floor(sigma*length(ind_deleted_edges));
[~,J] = maxk(X(:),k);
[i,j] = ind2sub([n n],J);
added = [i,j];

%% compute prediction accuracy (score)
score = ismember(added,G.Edges.EndNodes(ind_deleted_edges,:),'rows');
score = sum(score);

end




function v = nodangling_Tp(T,x,p)
    n = length(x);
    v1 = Tp(T,x,p);
    
    nonzer = reshape(double(ttv(T,{ones(n,1)},1)), [1,n^2]);
    d = ones(n^2,1) - nonzer';
    a = d'*reshape(meanp(x,x,p), [n^2,1]);
    v2 = (ones(n,1)./n)*a;
    
    v = v1 + v2;
end


function m = meanp(a,b,p)
    if p == 2
        u = abs(a).^2;
        v = abs(b).^2;
        m = sqrt((u + v') ./ 2);
    elseif p == 1
        m = (a./2) + (b./2)' ;
    elseif p == 0
        m = sqrt(a * b');
    else
        u = abs(a).^p;
        v = abs(b).^p;
        m = ((u + v') ./ 2).^(1/p);
    end
end


