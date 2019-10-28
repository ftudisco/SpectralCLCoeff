function Gmax = max_connected_subgraph(G)
    comp = conncomp(G);
    [~,val] = max(histc(comp,unique(comp)));
    nodes = 1 : numnodes(G);
    Gmax = subgraph(G, nodes(comp==val));
end