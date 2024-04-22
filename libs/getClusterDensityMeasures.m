function [ feat,featDescription ] = getClusterDensityMeasures(lympMembers,lympNodes,nonLympMembers,nonLympNodes)

[lympTot,lymAvg,lymMedian,lymStd,lymMode]=getDensityMeasures(lympMembers,lympNodes);
[nonLympTot,nonLymAvg,nonLymMedian,nonLymStd,nonLymMode]=getDensityMeasures(nonLympMembers,nonLympNodes);

feat=[lympTot,lymAvg,lymMedian,lymStd,lymMode,nonLympTot,nonLymAvg,...
    nonLymMedian,nonLymStd,nonLymMode];
featDescription={'lympTotDensity','lymAvgDensity','lymMedianDensity',...
    'lymStdDensity','lymModeDensity','nonLympTotDensity','nonLymAvgDensity',...
    'nonLymMedianDensity','nonLymStdDensity','nonLymModeDensity'};

end

function [totalDens,avgDens,medianDens,stdDens,modeDens]=getDensityMeasures(members,nodes)
numMembers=length(members);
dens=[];
for i = 1:numMembers
    member = members{i};
    col = nodes.centroid_c(member)';
    row = nodes.centroid_r(member)';
    numNodes=length(col);
    if  numNodes >  2
        [k,a] = convhull(col,row);
        dens=[dens;(numNodes/a)];
    end
end
totalDens=sum(dens);
totalDens(isnan(totalDens)) = 0;
avgDens=mean(dens);
avgDens(isnan(avgDens)) = 0;
medianDens=median(dens);
medianDens(isnan(medianDens)) = 0;
stdDens=std(dens);
stdDens(isnan(stdDens)) = 0;
modeDens=mode(dens);
modeDens(isnan(modeDens)) = 0;
end