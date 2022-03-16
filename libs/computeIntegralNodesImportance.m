function importance = computeIntegralNodesImportance(G,Nodes)

sCoord = Nodes(G.Edges.EndNodes(:,1),:);
tCoord = Nodes(G.Edges.EndNodes(:,2),:);

edgeD = diag(pdist2(sCoord,tCoord));

% edgeD = (edgeD+0.0001-min(edgeD))/(max(edgeD)-min(edgeD));
edgeD = edgeD./norm(edgeD,1);
weights = exp(-edgeD);
% figure
% [~,i] = sort(edgeD);
% plot(sort(edgeD),weights(i));%plot the weights distribution to 

WeightedDegree=[];
WeightedCloseness=[];
WeightedBetweenness=[];
WeightedEigenvector=[];

if ~isempty(weights)
    WeightedDegree = centrality(G,'degree','Importance',weights);%degree equal to the sum of weights of connected edges
    % plot(G,'XData',Nodes(:,1),'YData',Nodes(:,2),'EdgeLabel',weights,'NodeLabel',weightedDegree)    
    WeightedEigenvector = centrality(G,'eigenvector','Importance',weights);
    % plot(G,'XData',Nodes(:,1),'YData',Nodes(:,2),'EdgeLabel',weights,'NodeLabel',Weightedeigenvector)
end

if ~isempty(edgeD)
    WeightedCloseness = centrality(G,'closeness','Cost',edgeD);%closeness based on normalized edge distance, more central the point is, higher value
    % plot(G,'XData',Nodes(:,1),'YData',Nodes(:,2),'EdgeLabel',edgeD,'NodeLabel',WeightedCloseness)
    WeightedBetweenness = centrality(G,'betweenness','Cost',edgeD);
    % plot(G,'XData',Nodes(:,1),'YData',Nodes(:,2),'EdgeLabel',edgeD,'NodeLabel',WeightedBetweenness)
end

importance = struct('WeightedDegree',WeightedDegree,'WeightedCloseness',WeightedCloseness,...
    'WeightedBetweenness',WeightedBetweenness,'WeightedEigenvector',WeightedEigenvector);