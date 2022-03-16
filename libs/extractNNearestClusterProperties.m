function [ feat,featDescription ]  = extractNNearestClusterProperties( lympMembers,lympNodes,nonLympMembers,nonLympNodes, n  )
%EXTRACTTWONEARESTCLUSTERPROPERTIES Summary of this function goes here
%   Detailed explanation goes here

[lympClusterCent,~]=getClusterCentroids(lympMembers,lympNodes);
[nonLympClusterCent,~]=getClusterCentroids(nonLympMembers,nonLympNodes);

cent=[lympClusterCent;nonLympClusterCent];

%lbl=[ones(length(lympClusterCent),1);zeros(length(nonLympClusterCent),1)];

[numLympCent,~]=size(lympClusterCent);
[numNonLympCent,~]=size(nonLympClusterCent);
lbl=[ones(numLympCent,1);zeros(numNonLympCent,1)];

d=squareform(pdist(cent));
d(d==0)=Inf;

numClust=length(cent);
nearestNClusters=[];

if ~isempty(d)
    for i=1:numClust
        val=0;
        for j=1:n
            [~,ind]=min(d(i,:));
            d(i,ind)=Inf;
            val=val+lbl(ind);
        end    
        nearestNClusters=[nearestNClusters;val];
    end
end

avgNearest=mean(nearestNClusters);
avgNearest(isnan(avgNearest))=0;
medianNearest=median(nearestNClusters);
medianNearest(isnan(medianNearest))=0;
modeNearest=mode(nearestNClusters);
modeNearest(isnan(modeNearest))=0;
feat=[avgNearest,medianNearest,modeNearest];
featDescription={['avg_' n '_nearest'],['median_' n '_nearest'],['mode_' n '_nearest']};

end
