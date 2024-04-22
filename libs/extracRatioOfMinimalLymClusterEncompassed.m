function [feat,featDiscription] = extracRatioOfMinimalLymClusterEncompassed(LymPoints,NucleiPoints,...
    LymConcaveCenter,NucleiConcaveCenter)
feat = [];
featAreaRatio = [];
featDensityRatio = [];
if isempty(NucleiPoints) || length(NucleiPoints)<3
    feat = zeros(1,12);
else
    for i = 1:length(NucleiPoints)
        %imshow(lymNucleiOverlaid),hold on
        NuceliClusterPoints = NucleiPoints{i};
        if size(NuceliClusterPoints,1) > 2
            NuceliClusterCenter = NucleiConcaveCenter(i,:);
            
            [~,colNuc]=size(NuceliClusterCenter);
            [~,colLym]=size(LymConcaveCenter);
            
            [k,nucleiConvArea] = convhull(NuceliClusterPoints(:,1),NuceliClusterPoints(:,2));
            
            if colLym==colNuc
                distanceList = pdist2(NuceliClusterCenter,LymConcaveCenter);
                [~,Index] = sort(distanceList);
                %imshow(lymNucleiOverlaid),hold on
                %plot(NuceliClusterPoints(k,1),NuceliClusterPoints(k,2),'r','LineWidth',3);
                [LymTotalNum,LymTotalArea] = findMinAreaOfEncompassLymCluster(Index,LymConcaveCenter,NuceliClusterPoints,LymPoints);
            else
                LymTotalNum=0;
                LymTotalArea=0;
            end
            
            featAreaRatio = [featAreaRatio,nucleiConvArea/LymTotalArea];
            featDensityRatio = [featDensityRatio,(size(NuceliClusterPoints(:,1),1)/nucleiConvArea)/(LymTotalNum/LymTotalArea)];
            
        end
    end
    feat = [mean(featAreaRatio), std(featAreaRatio), var(featAreaRatio), min(featAreaRatio), max(featAreaRatio), max(featAreaRatio)-min(featAreaRatio),...
        mean(featDensityRatio), std(featDensityRatio), var(featDensityRatio), min(featDensityRatio), max(featDensityRatio), max(featDensityRatio)-min(featDensityRatio)];
end
feat(isnan(feat))=0;
feat(isinf(feat))=0;
if sum(feat) == 0
    feat = zeros(1,12);
end
featDiscription = {'avg of area ratio of nuclei cluster to encompassed lym clusters'...
    'std of area ratio of nuclei cluster to encompassed lym clusters'...
    'var of area ratio of nuclei cluster to encompassed lym clusters'...
    'min of area ratio of nuclei cluster to encompassed lym clusters'...
    'max of area ratio of nuclei cluster to encompassed lym clusters'...
    'range of area ratio of nuclei cluster to encompassed lym clusters'...
    'avg of density ratio of nuclei cluster to encompassed lym clusters'...
    'std of density ratio of nuclei cluster to encompassed lym clusters'...
    'var of density ratio of nuclei cluster to encompassed lym clusters'...
    'min of density ratio of nuclei cluster to encompassed lym clusters'...
    'max of density ratio of nuclei cluster to encompassed lym clusters'...
    'range of density ratio of nuclei cluster to encompassed lym clusters'};

end

function [LymTotalNum,LymTotalArea] = findMinAreaOfEncompassLymCluster(Index,LymConcaveCenter,NuceliClusterPoints,LymPoints)
EncompassVertices = [];
LymArea = [];
LymNumber = [];
LymTotalNum=0;
LymTotalArea=0;
for i = 1:length(Index)
    if size(LymPoints{Index(i)},1) >= 3
        
        
        
        EncompassVertices = [EncompassVertices;LymConcaveCenter(Index(i),:)];
        
        LymPoints{Index(i)};
        [k,area] = convhull(LymPoints{Index(i)}(:,1),LymPoints{Index(i)}(:,2));
        LymNumber =[LymNumber,size(LymPoints{Index(i)}(:,1),1)];
        LymArea = [LymArea, area];
        %plot(LymPoints{Index(i)}(k,1),LymPoints{Index(i)}(k,2),'b','LineWidth',3);
        
        
        if  size(EncompassVertices,1) == 1
            startpoint = LymConcaveCenter(Index(i),:);
            LymTotalArea = Inf;
            LymTotalNum = Inf;
        elseif size(EncompassVertices,1) > 2
            
            [k,~] = convhull(EncompassVertices(:,1),EncompassVertices(:,2));
            in = inpolygon(NuceliClusterPoints(:,1),NuceliClusterPoints(:,2),...
                [EncompassVertices(k,1);startpoint(1)],[EncompassVertices(k,2);startpoint(2)]);
            if all(in)
                %plot(NuceliClusterCenter(1),NuceliClusterCenter(2),'*r','LineWidth',3)
                %plot(EncompassVertices(k,1),EncompassVertices(k,2),'b','LineWidth',3);
                %close all
                LymTotalArea = sum(LymArea);
                LymTotalNum = sum(LymNumber);
                break
            elseif size(EncompassVertices,1) > 8
                LymTotalArea = Inf;
                LymTotalNum = Inf;
                break
            end
        end
        
    end
end

end