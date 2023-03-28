function [feats,feature_list,params] = cluster_graph_features(bounds, varargin)
% MATLAB Implementation of Sahir's ccg features
% Usage:
% cluster_graph_features_networkboost(bounds)
% cluster_graph_features_networkboost(bounds, edges)
% cluster_graph_features_networkboost(bounds, alpha, r)
%
% updated (Ajay, 2014) - additional input options allow CCG to be constructed on-the-fly
% updated (Ajay, 2014) - reduced complexity from O(n^2) to 0(n)
% updated - faster implementation of connected components
% George Lee (2013)

% Allows for 3 different input argument structures 
params = [];
if nargin == 1
    [~, ~, ~, ~, edges, params] = construct_ccgs(bounds, 0.45, 0.21);
elseif nargin == 2
    edges = varargin{1};
elseif nargin == 3
    [~, ~, ~, ~, edges, params] = construct_ccgs(bounds, varargin{1}, varargin{2});
else 
    error('Input args must be: (bounds), (bounds, edges), or (bounds, alpha, r)');
end

%1) Number of Nodes
N = length(edges);
feats(1) = N;
feature_list(1) = {'Number of Nodes'};


%2) Number of edges
E = nnz(edges);
feats(2) = E;
feature_list(2) = {'Number of Edges'};


%3) Average Degree
feats(3) = E/N;
feature_list(3) = {'Average Degree'};


%%% Eccentricity calculation
% generate distance matrix
X = [bounds(:).centroid_r; bounds(:).centroid_c];
D = squareform(pdist(X','euclidean')) + eye(length(X));

% create sparse distance-weighted edge matrix
edges = triu(edges); % force edges to be upper triangular
sym_edges = sparse(edges | edges'); % full symmetric matrix
weighted = sym_edges .* D; 

% pre-calculate shortest paths (8/24/2013)
for i = 1:N
    distance = graphshortestpath(weighted,i);
    pathlengths{i} = nonzeros(distance(isfinite(distance)));
end
pathlengths_nonempty = pathlengths;
pathlengths_nonempty(cellfun(@isempty, pathlengths_nonempty, 'UniformOutput', true)) = [];

% all non-zero pathlengths
eccentricity = cellfun(@max, pathlengths_nonempty, 'UniformOutput', true);

%4) Average eccentricity
feats(4) = sum(eccentricity)/N;
feature_list(4) = {'Average Eccentricity'};

%5) Diameter
try
    diameter = max(eccentricity);
    feats(5) = diameter;
    feature_list(5) = {'Diameter'};
catch
    feats(5)=NaN;
    feature_list(5) = {'Diameter'};
end

%6) Radius
try
    radius = min(eccentricity);
    feats(6) = radius;
    feature_list(6) = {'Radius'};
catch
    feats(6)=NaN;
    feature_list(6) = {'Radius'};
end


% eccentricity for largest 90% of path lengths

%%% CHOOSE ONE OF THE FOLLOWING %%%
%%%%% George's definition of 90% %%%%%
pathlengths_nonempty90 = cellfun(@(x) x(1:round(0.9*length(x))), pathlengths_nonempty, 'UniformOutput', false);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Ajay's definition of 90% %%%%%
% pathlengths_nonempty90 = cellfun(@(x) x(x >= prctile(x,10)), pathlengths_nonempty, 'UniformOutput', false);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pathlengths_nonempty90(cellfun(@isempty, pathlengths_nonempty90, 'UniformOutput', true)) = [];
eccentricity90 = cellfun(@max, pathlengths_nonempty90, 'UniformOutput', true);

%7) Average Eccentricity 90%
try
    feats(7) = sum(eccentricity90)/N;   
    feature_list(7) = {'Average Eccentricity 90\%'};
catch
    feats(7)=NaN;
    feature_list(7) = {'Average Eccentricity 90\%'};
end


%8) Diameter 90%
try
    diameter90 = max(eccentricity90);
    feats(8) = diameter90;
    feature_list(8) = {'Diameter 90\%'};
catch
    feats(8) = NaN;
    feature_list(8) = {'Diameter 90\%'};
end

%9) Radius 90%
try
    radius90 = min(eccentricity90);
    feats(9) = radius90;
    feature_list(9) = {'Radius 90\%'};
catch
    feats(9) = NaN;
    feature_list(9) = {'Radius 90\%'};
end
    
%10) Average Path Length
try

    feats(10) = sum(cellfun(@sum, pathlengths, 'UniformOutput', true)) / sum(cellfun(@length, pathlengths, 'UniformOutput', true));
    feature_list(10) = {'Average Path Length'};
catch
    feats(10) = sum(cellfun(@sum, pathlengths, 'UniformOutput', true)) / sum(cellfun(@length, pathlengths, 'UniformOutput', true));
    feature_list(10) = {'Average Path Length'};
end


%% clustering coefficients
[~,network] = graphconncomp(sym_edges);

for n=1:N
    nodes = find(network==n);
    En(nodes) = sum(sum(edges(nodes, nodes)));
    kn(nodes) = length(nodes);
end

%11) Clustering Coefficient C
% ratio beteween A: the number of edges between neighbors of node n and B:
% the number of possible edges between the neighbors of node n

%Cn(n) = 2*En(n) / ( kn(n) * (kn(n)-1) );
Cn = 2*En ./ ( kn .* (kn-1) );
Cn( isnan(Cn) ) = 0; % account for divide by zero
feats(11) = sum(Cn)/N;
feature_list(11) = {'Clustering Coefficient C'};


%12) Clustering Coefficient D
%Dn(n) = 2*(kn(n) + En(n)) / ( kn(n)* (kn(n)+1) );
Dn = 2*(kn + En) ./ ( kn .* (kn+1) );
Dn( isnan(Dn) ) = 0;
feats(12) = sum(Dn)/N;
feature_list(12) = {'Clustering Coefficient D'};


%13) Clustering Coefficient E
% count isolated nodes
iso_nodes = sum(kn == 1);
feats(13) = sum( Cn(kn > 1) ) / (N - iso_nodes);
feature_list(13) = {'Clustering Coefficient E'};

%14) Number of connected components
feats(14) = length(kn(kn>1));
feature_list(14) = {'Number of connected components'};


%15) Giant connected component ratio
feats(15) = max(kn) / N;
feature_list(15) = {'giant connected component ratio'};

%16) Average Connected Component Size
feats(16) = mean(kn(kn>1));
feature_list(16) = {'average connected component size'};


%17 and 18) Number / Percentage of Isolated Nodes
feats(17) = iso_nodes;
feature_list(17) = {'number isolated nodes'};

feats(18) = iso_nodes/N;
feature_list(18) = {'percentage isolated nodes'};


%19 and 20) Number / Percentage of End points
feats(19) = sum(kn==2);
feature_list(19) = {'number end nodes'};
feats(20) = sum(kn==2)/N;
feature_list(20) = {'percentage end nodes'};


% 21 and 22) Number / Percentage of Central points
feats(21) = sum(eccentricity == radius);
feature_list(21) = {'number central nodes'};
feats(22) = sum(eccentricity == radius) / N;
feature_list(22) = {'percentage central nodes'};

% 23 - 26) Edge length statistics
edge_lengths = weighted(:);
edge_lengths(edge_lengths==0) = []; % remove zero edge lengths

feats(23) = sum(edge_lengths)/length(edge_lengths); % mean edge-length
feature_list(23) = {'mean edge length'};
feats(24) = std(edge_lengths); % standard deviation
feature_list(24) = {'standard deviation edge length'};
feats(25) = skewness(edge_lengths); % skewness
feature_list(25) = {'skewness edge length'};
feats(26) = kurtosis(edge_lengths); % kurtosis
feature_list(26) = {'kurtosis edge length'};