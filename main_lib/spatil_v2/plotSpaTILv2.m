function plotSpaTILv2(args)
%PLOTSPATILV2 Summary of this function goes here
%   Detailed explanation goes here


I=args.I;
maskNucFile=args.maskNucFile;
maskEpiFile=args.maskEpiFile;
lympDataFile=args.lympDataFile;
plot=args.plot;


if ~isfield(args,'alpha')
    alpha=.45;
else
    alpha=args.alpha;
end

if ~isfield(args,'r')
    r=.185;
else
    r=args.r;
end



%% Showing original image

if strcmp(plot,'image')
    imshow(I); % Shows the tile
    return;
end

%% Nuclei mask

MN=imread(maskNucFile);
if size(MN,3)>1
    MN(:,:,1:2)=[];
end
MN=imbinarize(MN,graythresh(MN));

if strcmp(plot,'nuclei_alpha')
    
    %imshow(imoverlay(I,MN,[0 1 0])); % Overlays nuclei (no alpha)
    imshow(I); alphamask(logical(MN),[0 0 1],.1); % Overlays nuclei (with alpha)
    return;
end

%% Epithelium mask

if strcmp(plot,'epithelium')
    
    ME=imread(maskEpiFile);
    if size(ME,3)>1
        ME=rgb2gray(ME);
    end
    %ME=imresize(ME,4);
    
    imshow(I); alphamask(logical(ME),[0 1 1],.3); % Overlays epithelium (with alpha)
    return;
end

%% Reading lymphocyte-related data
data=load(lympDataFile);
classes=getEpiStrLympClasses(data.isLymphocyte,~data.epiNuc);
numCent=length(data.nucCentroids);

if strcmp(plot,'nuclei_borders')
    drawNucContoursByClass(MN,I,data.nucCentroids,ones(numCent,1)*4,5); % Plots nuclei borders using a single color
    return;
end

if strcmp(plot,'nuclei_dots')
    drawNucleiCentroidsByClass(I,data.nucCentroids,ones(numCent,1)*3,8); % Plots nuclei as dots using a single color
    return;
end

if strcmp(plot,'nuclei_dots_multi')
    drawNucleiCentroidsByClass(I,data.nucCentroids,classes,10); % Plots nuclei as dots with diff. color per class.
    return;
end

if strcmp(plot,'nuclei_boders_multi')
    drawNucContoursByClass(MN,I,data.nucCentroids,classes,5); % Plots nuclei borders with diff. color per class.
    return;
end

%% Cell cluster graphs

coords= {
    data.nucCentroids(~data.isLymphocyte & data.epiNuc,:),...
    data.nucCentroids(data.isLymphocyte & ~data.epiNuc,:),...
    data.nucCentroids(data.isLymphocyte & data.epiNuc,:),...
    data.nucCentroids(~data.isLymphocyte & ~data.epiNuc,:),...
    };

numGroups=length(coords);

M=cell(numGroups,1);
for i=1:numGroups
    [~,~,~,~,groupMatrix] = construct_nodesCluster(struct('centroid_r',coords{i}(:,2)','centroid_c',coords{i}(:,1)'), alpha, r);
    M{i}=groupMatrix;
end

if strcmp(plot,'graphs')
    %drawGraph(I,coords,M);
    drawGraph(I,coords,M,10,100,.8,.8,true);
    return;
end

if strcmp(plot,'clusters')
    
    args.fillAlpha=.3;
    args.lineWidth=10;
    args.lineAlpha=.7;
    args.useBlackLine=true;
    
    plotClusters(I,coords,M,args);
    %plotClusters(I,coords,M,true,0.5,5,0,false);
    return;
end

end

