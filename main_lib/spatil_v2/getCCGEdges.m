function edges = getCCGEdges(nodes,alpha,r)
%GETCCGEDGES Summary of this function goes here
%   Detailed explanation goes here
X = [nodes(:).centroid_r; nodes(:).centroid_c];

D = squareform(pdist(X','euclidean'));
P = D.^-alpha;

edges = triu(P>r,1);
