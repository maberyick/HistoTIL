function [feat,featDiscription] = extractLocalMinimalLymClusterEncompassed(LymPoints,NucleiPoints,...
    LymConcaveCenter,NucleiConcaveCenter)
feat = [];
if isempty(NucleiPoints) || length(NucleiPoints)<3
    feat = zeros(1,6);
else
    %for i=4%i = 1:length(NucleiPoints)%[13] [4]
    %imshow(lymNucleiOverlaid),hold on
    %%
    if length(NucleiPoints)>3
        i=4;
        NuceliClusterPoints = NucleiPoints{i};
        if size(NuceliClusterPoints,1) > 2
            NuceliClusterCenter = NucleiConcaveCenter(i,:);
            [~,colNuc]=size(NuceliClusterCenter);
            [~,colLym]=size(LymConcaveCenter);
            if colLym==colNuc
                distanceList = pdist2(NuceliClusterCenter,LymConcaveCenter);
                [~,Index] = sort(distanceList);
                %imshow(lymNucleiOverlaid),hold on
                %[k,nucleiConvArea] = convhull(NuceliClusterPoints(:,1),NuceliClusterPoints(:,2));
                %plot(NuceliClusterPoints(k,1),NuceliClusterPoints(k,2),'y','LineWidth',3);
                %h = fill(NuceliClusterPoints(k,1),NuceliClusterPoints(k,2),[1 1 0]);
                %set(h,'facealpha',0.5)
                [LymNumber] = findMinNumberOfEncompassLym(NuceliClusterCenter,Index,LymConcaveCenter,NuceliClusterPoints,LymPoints);
            else
                LymNumber=0;
            end
            feat = [feat,1/(LymNumber-2)];
            feat(isnan(feat))=0;
            if isempty(feat)
                feat = 0;
            end
        end
        %%
    end
    %end
    %close all
    feat(isnan(feat))=0;
    if isempty(feat)
        feat = 0;
    end
    featV = feat;
    feat = [mean(featV), std(featV), var(featV), min(featV), max(featV), max(featV)-min(featV)];
end
featDiscription = {'avg of reciprocal of number of least lym cluster to encompass nuclei cluster'...
    'std of reciprocal of number of least lym cluster to encompass nuclei cluster'...
    'var of reciprocal of number of least lym cluster to encompass nuclei cluster'...
    'min of reciprocal of number of least lym cluster to encompass nuclei cluster'...
    'max of reciprocal of number of least lym cluster to encompass nuclei cluster'...
    'range of reciprocal of number of least lym cluster to encompass nuclei cluster'};

end

function [LymNumber] = findMinNumberOfEncompassLym(NuceliClusterCenter,Index,LymConcaveCenter,NuceliClusterPoints,LymPoints)
EncompassVertices = [];
LymNumber = Inf;
for j = 1:length(Index)
    if size(LymPoints{Index(j)},1) >= 3
        %[k,~] = convhull(LymPoints{Index(j)}(:,1),LymPoints{Index(j)}(:,2));
        %plot(LymPoints{Index(j)}(k,1),LymPoints{Index(j)}(k,2),'b','LineWidth',1);
        %h = fill(LymPoints{Index(j)}(k,1),LymPoints{Index(j)}(k,2),[0 1 0]);
        %set(h,'facealpha',0.3)
        EncompassVertices = [EncompassVertices;LymConcaveCenter(Index(j),:)];
        if  size(EncompassVertices,1) == 1
            startpoint = LymConcaveCenter(Index(j),:);
            LymNumber = Inf;
        elseif size(EncompassVertices,1) > 2
            
            [k,~] = convhull(EncompassVertices(:,1),EncompassVertices(:,2));
            in = inpolygon(NuceliClusterPoints(:,1),NuceliClusterPoints(:,2),...
                [EncompassVertices(k,1);startpoint(1)],[EncompassVertices(k,2);startpoint(2)]);
            if all(in)
                %plot([EncompassVertices(k,1);startpoint(1)],[EncompassVertices(k,2);startpoint(2)])
                %plot(NuceliClusterCenter(1),NuceliClusterCenter(2),'*y','LineWidth',3)
                %plot(EncompassVertices(k,1),EncompassVertices(k,2),'b','LineWidth',3);
                %                 close all
                LymNumber = size(EncompassVertices,1);
                break
            elseif size(EncompassVertices,1) > 6
                LymNumber = Inf;
                break
            end
        end
    end
end

end