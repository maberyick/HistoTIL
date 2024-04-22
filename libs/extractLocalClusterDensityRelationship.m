function [feat,featDiscription] = extractLocalClusterDensityRelationship(LymPoints,NucleiPoints,...
    LymConcaveCenter,NucleiConcaveCenter)
featV = [];
if isempty(NucleiPoints) || length(NucleiPoints)<3
    feat = zeros(1,24);
else
    %%
    for i = 1:length(NucleiPoints)
        %figure,imshow(lymNucleiOverlaid),hold on
        NuceliClusterPoints = NucleiPoints{i};
        if size(NuceliClusterPoints,1) > 2
            NuceliClusterCenter = NucleiConcaveCenter(i,:);
            [~,colNuc]=size(NuceliClusterCenter);
            [~,colLym]=size(LymConcaveCenter);
            if colLym==colNuc
                distanceList = pdist2(NuceliClusterCenter,LymConcaveCenter);
                [~,Index] = sort(distanceList);
                [lymDensity] = findNearest3LymCluster(Index,LymPoints);
                %                 LymPonitsIn = LymPoints{LymIndex};
                [k,nucleiConvArea] = convhull(NuceliClusterPoints(:,1),NuceliClusterPoints(:,2));
                %plot(NuceliClusterPoints(k,1),NuceliClusterPoints(k,2),'r','LineWidth',3);
                nucleiDensity = size(NuceliClusterPoints,1)/nucleiConvArea;
                featEach = [nucleiDensity./lymDensity,nucleiDensity/mean(lymDensity)];
                featEach(isnan(featEach))=0;
            else
                featEach=0;
            end
            featV = [featV;featEach];
            %close all
        end
    end
    %%
    %mean
    if  isempty(featV)
        feat = zeros(1,24);
    elseif featV==0
        feat = zeros(1,24);
    elseif length(featV)<24
        feat = zeros(1,24);
    else
        avgfeat = mean(featV);
        %std
        stdfeat = std(featV);
        %var
        varfeat = var(featV);
        %min
        minfeat = min(featV);
        %max
        maxfeat = max(featV);
        %range
        rangefeat = maxfeat-minfeat;
        %feat
        feat = [avgfeat,stdfeat,varfeat,minfeat,maxfeat,rangefeat];
    end
end
%featDiscription
featDiscription = {'avg of ratio of the nuceli cluster density to its closest lym cluster cluster',...
    'avg of ratio of the nuceli cluster density to its second closest lym cluster density',...
    'avg of ratio of the nuceli cluster density to its third closest lym cluster density',...
    'avg of ratio of the nuceli cluster density to its first three closest lym cluster mean density',...
    'std of ratio of the nuceli cluster density to its closest lym cluster cluster',...
    'std of ratio of the nuceli cluster density to its second closest lym cluster density',...
    'std of ratio of the nuceli cluster density to its third closest lym cluster density',...
    'std of ratio of the nuceli cluster density to its first three closest lym cluster mean density',...
    'var of ratio of the nuceli cluster density to its closest lym cluster cluster',...
    'var of ratio of the nuceli cluster density to its second closest lym cluster density',...
    'var of ratio of the nuceli cluster density to its third closest lym cluster density',...
    'var of ratio of the nuceli cluster density to its first three closest lym cluster mean density',...
    'min of ratio of the nuceli cluster density to its closest lym cluster cluster',...
    'min of ratio of the nuceli cluster density to its second closest lym cluster density',...
    'min of ratio of the nuceli cluster density to its third closest lym cluster density',...
    'min of ratio of the nuceli cluster density to its first three closest lym cluster mean density',...
    'max of ratio of the nuceli cluster density to its closest lym cluster cluster',...
    'max of ratio of the nuceli cluster density to its second closest lym cluster density',...
    'max of ratio of the nuceli cluster density to its third closest lym cluster density',...
    'max of ratio of the nuceli cluster density to its first three closest lym cluster mean density',...
    'range of ratio of the nuceli cluster density to its closest lym cluster cluster',...
    'range of ratio of the nuceli cluster density to its second closest lym cluster density',...
    'range of ratio of the nuceli cluster density to its third closest lym cluster density',...
    'range of ratio of the nuceli cluster density to its first three closest lym cluster mean density'};
end

function [lymDensity] = findNearest3LymCluster(Index,LymPoints)
LymIndex =[];
lymDensity = [];
for i = 1:length(Index)
    if size(LymPoints{Index(i)},1) >= 3
        LymIndex = [LymIndex,Index(i)];
        [k,area] = convhull(LymPoints{Index(i)}(:,1),LymPoints{Index(i)}(:,2));
        plot(LymPoints{Index(i)}(k,1),LymPoints{Index(i)}(k,2),'b','LineWidth',3);
        lymDensity = [lymDensity,size(LymPoints{Index(i)},1)/area];
    end
    if length(LymIndex) == 3
        break;
    end
end
if length(lymDensity) == 2
    lymDensity = [lymDensity,0];
elseif length(lymDensity) == 1
    lymDensity = [lymDensity,0,0];
end
end