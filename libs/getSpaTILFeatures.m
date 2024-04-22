function [features] = getSpaTILFeatures(lympCentroids_cur,nonLympCentroids_cur)
%GETSPATILFEATURES Gets a set of spatial quantitative features (SpaTIL) 
% relating to the spatial architecture of TILs, co-localization of TILs and 
% cancer nuclei, and density variation of TIL clusters.
alpha=.49;
r=.185;

lympNodes = struct('centroid_r',lympCentroids_cur(:,2)','centroid_c',lympCentroids_cur(:,1)');
nonLympNodes = struct('centroid_r',nonLympCentroids_cur(:,2)','centroid_c',nonLympCentroids_cur(:,1)');

[~,~,~,~,lympMatrix,~,~,~] = construct_nodesCluster(lympNodes, alpha, r);
[~,~,~,~,nonLympMatrix,~,~,~] = construct_nodesCluster(nonLympNodes, alpha, r);

[~,~,lympMembers] = networkComponents(lympMatrix);
[~,~,nonLympMembers] = networkComponents(nonLympMatrix);

[ feat5,~ ] = getClusterAreaMeasures(lympMembers,lympNodes,nonLympMembers,nonLympNodes);
[ feat6,~] = getClusterDensityMeasures(lympMembers,lympNodes,nonLympMembers,nonLympNodes);
[ feat7,~] = getIntersectedArea( lympMembers,lympNodes,nonLympMembers,nonLympNodes );

[ feat8_1,~ ]  = extractNNearestClusterProperties(lympMembers,lympNodes,nonLympMembers,nonLympNodes,1);
[ feat8_2,~ ]  = extractNNearestClusterProperties( lympMembers,lympNodes,nonLympMembers,nonLympNodes,2);
[ feat8_3,~ ]  = extractNNearestClusterProperties( lympMembers,lympNodes,nonLympMembers,nonLympNodes,3);
[ feat8_4,~ ]  = extractNNearestClusterProperties( lympMembers,lympNodes,nonLympMembers,nonLympNodes,4);
[ feat8_5,~ ]  = extractNNearestClusterProperties( lympMembers,lympNodes,nonLympMembers,nonLympNodes,5);
feat8 = [feat8_1,feat8_2,feat8_3,feat8_4,feat8_5];
%featDescription8 = [featDescription81,featDescription82,featDescription83,featDescription84,featDescription85];

Lymnodes = struct('centroid_r',lympCentroids_cur(:,2)','centroid_c',lympCentroids_cur(:,1)');
Nucleinodes = struct('centroid_r',nonLympCentroids_cur(:,2)','centroid_c',nonLympCentroids_cur(:,1)');

[~,~,~,~,s_edges] = construct_nodesCluster(Lymnodes, 0.49, 0.185);
[~,~,~,~,s_edgesNuclei] = construct_nodesCluster(Nucleinodes, 0.49, 0.185);

s_edges = s_edges|s_edges';
LymConcaveCenter = [];
[~,~,Lymmembers] = networkComponents(s_edges);
LymPoints = {};
numLymMembers=length(Lymmembers);
parfor j = 1:numLymMembers
    a = Lymmembers{j};
    LymConcaveC = Lymnodes.centroid_c(a)';
    LymConcaveR = Lymnodes.centroid_r(a)';
    LymPoints(j) = {[LymConcaveC,LymConcaveR]};
    ConcaveCenter = [mean(LymConcaveC),mean(LymConcaveR)];
    LymConcaveCenter = [LymConcaveCenter;ConcaveCenter];
end

s_edgesNuclei = s_edgesNuclei|s_edgesNuclei';
NucleiConcaveCenter = [];
[~,~,Nucelimembers] = networkComponents(s_edgesNuclei);
NucleiPoints = {};
numNucleiMembers=length(Nucelimembers);
parfor j = 1:numNucleiMembers
    a = Nucelimembers{j};
    NucleiConcaveC = Nucleinodes.centroid_c(a)';
    NucleiConcaveR = Nucleinodes.centroid_r(a)';
    NucleiPoints(j) = {[NucleiConcaveC,NucleiConcaveR]};
    ConcaveCenter = [mean(NucleiConcaveC),mean(NucleiConcaveR)];
    NucleiConcaveCenter = [NucleiConcaveCenter;ConcaveCenter];
end

[feat1,~] = extractLocalClusterDensityRelationship(LymPoints,NucleiPoints,LymConcaveCenter,NucleiConcaveCenter);
[feat2,~] = extractLocalMinimalLymClusterEncompassed(LymPoints,NucleiPoints,LymConcaveCenter,NucleiConcaveCenter);
[feat3,~] = extracRatioOfMinimalLymClusterEncompassed(LymPoints,NucleiPoints,LymConcaveCenter,NucleiConcaveCenter);
[NucleiedgesConnection,NuceliClusterCenterKept,LymedgesConnection,LymClusterCenterKept] = extractClusterGraph(LymPoints,NucleiPoints,LymConcaveCenter,NucleiConcaveCenter);

a = logical(NucleiedgesConnection);
a = tril(a)'|a|triu(a)';
NucleiedgesConnection = a;
NucleiG = graph(NucleiedgesConnection);
a = logical(LymedgesConnection);
a = tril(a)'|a|triu(a)';
LymedgesConnection = a;
LymG = graph(LymedgesConnection);
NucleiImportance = computeIntegralNodesImportance(NucleiG,NuceliClusterCenterKept);
LymImportance = computeIntegralNodesImportance(LymG,LymClusterCenterKept);
[feat4,~] = findIndexOfNearestLym(NucleiImportance,NuceliClusterCenterKept,LymImportance,LymClusterCenterKept);

features = [feat1,feat2,feat3,feat4,feat5,feat6,feat7,feat8];
close all;
%featureNames = [featDescription1,featDescription2,featDescription3,featDescription4,featDescription5,featDescription6,featDescription7,featDescription8];
end