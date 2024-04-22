function [NucleiedgesConnection,NuceliClusterCenterKept,LymedgesConnection,LymClusterCenterKept] = extractClusterGraph(LymPoints,NucleiPoints,LymConcaveCenter,NucleiConcaveCenter)

% close all,
% hold on
% imshow(lymNucleiOverlaid),hold on

NuceliClusterPointsKept = {};
NuceliClusterCenterKept = [];
VNucleiConcaveArea = [];
VNucleiNumberCluster = [];
for i = 1:length(NucleiPoints)
    NuceliClusterPoints = NucleiPoints{i};
    nucleiConvArea = 0;
    if size(NuceliClusterPoints,1) > 2
        %i
%         NuceliClusterCenter = NucleiConcaveCenter(i,:);
%         distanceList = pdist2(NuceliClusterCenter,LymConcaveCenter);
%         [~,Index] = sort(distanceList);
        %         imshow(lymNucleiOverlaid),hold on
        [k,nucleiConvArea] = convhull(NuceliClusterPoints(:,1),NuceliClusterPoints(:,2));
%                 plot(NuceliClusterPoints(k,1),NuceliClusterPoints(k,2),'r','LineWidth',3);
        
%         plot(NuceliClusterPoints(k,1),NuceliClusterPoints(k,2),'r','LineWidth',1);
        NuceliClusterPointsKept = [NuceliClusterPointsKept,[NuceliClusterPoints(k,1),NuceliClusterPoints(k,2)]];
        NuceliClusterCenterKept = [NuceliClusterCenterKept;NucleiConcaveCenter(i,:)];
        VNucleiConcaveArea = [VNucleiConcaveArea, nucleiConvArea];
        VNucleiNumberCluster = [VNucleiNumberCluster, size(NuceliClusterPoints,1)];
    end    
end

MDNucleiConcaveCenter = pdist(NuceliClusterCenterKept,'euclidean');
MDNucleiConcaveCenter = squareform(MDNucleiConcaveCenter);
MNucleiArea = [];
MNucleiDensity = [];
for i = 1:length(MDNucleiConcaveCenter)-1%only need to compute the distance of second last with the last, no needs to iterate to last, so length()-1
%     count = 0;
    for j = i+1:length(MDNucleiConcaveCenter) %do not need to compute self-to-self, so starting from the next one
        MNucleiArea(i,j) = VNucleiConcaveArea(i) + VNucleiConcaveArea(j);
        MNucleiDensity(i,j) = VNucleiNumberCluster(i)/VNucleiConcaveArea(i) + VNucleiNumberCluster(j)/VNucleiConcaveArea(j);
                
    end
end

MNucleiArea=0;
MNucleiDensity=0;
MDNucleiConcaveCenter=0;

if ~isempty(MNucleiArea)
    MNucleiArea = (MNucleiArea - min(MNucleiArea(:)))/(max(MNucleiArea(:))-min(MNucleiArea(:)));MNucleiArea = [MNucleiArea;zeros(1,size(MNucleiArea,2))];
end

if ~isempty(MNucleiDensity)
    MNucleiDensity = (MNucleiDensity - min(MNucleiDensity(:)))/(max(MNucleiDensity(:))-min(MNucleiDensity(:)));MNucleiDensity = [MNucleiDensity;zeros(1,size(MNucleiDensity,2))];    
end

if ~isempty(MDNucleiConcaveCenter)
    MDNucleiConcaveCenter = (MDNucleiConcaveCenter - min(MDNucleiConcaveCenter(:)))/(max(MDNucleiConcaveCenter(:))-min(MDNucleiConcaveCenter(:)));
end

MNucleiComb = 0.1.*MNucleiArea + 0.2.*MNucleiDensity + 0.7.*MDNucleiConcaveCenter;
NucleiedgesConnection=[];
if ~isnan(MNucleiComb)
    MNucleiComb = triu(MNucleiComb)';
    [ST,pred] = graphminspantree(sparse(MNucleiComb));
    NucleiedgesConnection = full(ST);
end
% gplot(NucleiedgesConnection,[NuceliClusterCenterKept(:,1),NuceliClusterCenterKept(:,2)])


% hold on
% imshow(lymNucleiOverlaid),hold on
LymClusterPointsKept = {};
LymClusterCenterKept = [];
VLymConcaveArea = [];
VLymNumberCluster = [];
for i = 1:length(LymPoints)
    LymClusterPoints = LymPoints{i};
    LymConvArea = 0;
    if size(LymClusterPoints,1) > 2
        %i
%         LymClusterCenter = LymConcaveCenter(i,:);
%         distanceList = pdist2(LymClusterCenter,NucleiConcaveCenter);
%         [~,Index] = sort(distanceList);
        %         imshow(lymNucleiOverlaid),hold on
        [k,LymConvArea] = convhull(LymClusterPoints(:,1),LymClusterPoints(:,2));
        %         plot(NuceliClusterPoints(k,1),NuceliClusterPoints(k,2),'r','LineWidth',3);
        
%         plot(LymClusterPoints(k,1),LymClusterPoints(k,2),'b','LineWidth',1);
        LymClusterPointsKept = [LymClusterPointsKept,[LymClusterPoints(k,1),LymClusterPoints(k,2)]];
        LymClusterCenterKept = [LymClusterCenterKept;LymConcaveCenter(i,:)];
        VLymConcaveArea = [VLymConcaveArea, LymConvArea];
        VLymNumberCluster = [VLymNumberCluster, size(LymClusterPoints,1)];
    end    
end

MDLymConcaveCenter = pdist(LymClusterCenterKept,'euclidean');
MDLymConcaveCenter = squareform(MDLymConcaveCenter);
MLymArea = [];
MLymDensity = [];

if ~isempty(MDLymConcaveCenter)
    for i = 1:length(MDLymConcaveCenter)-1%only need to compute the distance of second last with the last, no needs to iterate to last, so length()-1
    %     count = 0;
        for j = i+1:length(MDLymConcaveCenter) %do not need to compute self-to-self, so starting from the next one
            MLymArea(i,j) = VLymConcaveArea(i) + VLymConcaveArea(j);
            MLymDensity(i,j) = VLymNumberCluster(i)/VLymConcaveArea(i) + VLymNumberCluster(j)/VLymConcaveArea(j);

        end
    end
end

if ~isempty(MLymArea)
    MLymArea = (MLymArea - min(MLymArea(:)))/(max(MLymArea(:))-min(MLymArea(:)));MLymArea = [MLymArea;zeros(1,size(MLymArea,2))];
end

if ~isempty(MLymDensity)
    MLymDensity = (MLymDensity - min(MLymDensity(:)))/(max(MLymDensity(:))-min(MLymDensity(:)));MLymDensity = [MLymDensity;zeros(1,size(MLymDensity,2))];
end

if ~isempty(MDLymConcaveCenter)
    MDLymConcaveCenter = (MDLymConcaveCenter - min(MDLymConcaveCenter(:)))/(max(MDLymConcaveCenter(:))-min(MDLymConcaveCenter(:)));
end


MLymComb = 0.1.*MLymArea + 0.2.*MLymDensity + 0.7.*MDLymConcaveCenter;
LymedgesConnection=[];
if ~isnan(MNucleiComb)
    
    MLymComb = triu(MLymComb)';
    [ST,pred] = graphminspantree(sparse(MLymComb));
    LymedgesConnection = full(ST);
end
% gplot(LymedgesConnection,[LymClusterCenterKept(:,1),LymClusterCenterKept(:,2)])

%[i,j] = find(LymedgesConnection);
%[~, p] = sort(max(i,j));
%i = i(p);
%j = j(p);
%XLym = [ LymClusterCenterKept(i,1) LymClusterCenterKept(j,1)]';
%YLym = [ LymClusterCenterKept(i,2) LymClusterCenterKept(j,2)]';
% figure,plot(XLym,YLym)

%[i,j] = find(NucleiedgesConnection);
%[~, p] = sort(max(i,j));
%i = i(p);
%j = j(p);
%XNuclei = [ NuceliClusterCenterKept(i,1) NuceliClusterCenterKept(j,1)]';
%YNuclei = [ NuceliClusterCenterKept(i,2) NuceliClusterCenterKept(j,2)]';

% figure,plot(XNuclei,YNuclei)

a = logical(LymedgesConnection);
a = tril(a)'|a|triu(a)';
LymedgesConnection = a;

if size(LymClusterCenterKept,1) > size(NuceliClusterCenterKept,1)
    largegraphAdj = LymedgesConnection;
    largegraphNodes = LymClusterCenterKept;
    numNeedsCut = size(LymClusterCenterKept,1) - size(NuceliClusterCenterKept,1);
elseif size(LymClusterCenterKept,1) == size(NuceliClusterCenterKept,1)%needing put pruningLargeTree to if at the end to control equal situation
    return
else
    largegraphAdj = NucleiedgesConnection;
    largegraphNodes = NuceliClusterCenterKept;
    numNeedsCut = size(NuceliClusterCenterKept,1) - size(LymClusterCenterKept,1);
end



% [newAdj,newPoints] = pruningLargeTree(largegraphAdj,largegraphNodes,numNeedsCut);
% [K, S, best_cost] = shape_matching(newPoints, NuceliClusterCenterKept, 'aco','shape_context','chisquare');
% S=graph_similarity(logical(NucleiedgesConnection),logical(NucleiedgesConnection));

% plot(LymConcaveCenter(:,1),LymConcaveCenter(:,2))

feat = [];
featDiscription = [];
end
