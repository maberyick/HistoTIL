clear;clc;

%% Reading data 

%%-- Option 1: from WSI
 
 %imgName='20150528_163928';
 imgName='20150528_162253';
 %imgName='20150528_173858';
 
 imgExt='.tiff';
 imgFolder='Z:\data\JHU_WU_Oropharyngeal_WSI_ITH\';
 featFolder='D:\German\Data\Oroph_WSI_JHU\';
 
 tileSize=2048;
 %x=50760;y=89400; % example
 x=18433; y=47105; % good prognosis
 %x=14337; y=16385; % poor prognosis
 
imgFile=[imgFolder imgName imgExt];
info=imfinfo(imgFile);
[~,ind] = max(cat(1,info.Height));
 
 I=imread(imgFile,'PixelRegion',{[x x+tileSize-1],[y y+tileSize-1]},ind);
 
 maskNucFile=[featFolder 'masks\MN_seg\' imgName '\' imgName '_x' num2str(x) '_y' num2str(y) '_result.png'];
 maskEpiFile=[featFolder 'masks\epi_seg\' imgName '\' imgName '_x' num2str(x) '_y' num2str(y) '_result.png'];
 lympDataFile=[featFolder 'features\TIL_detection\' imgName '\' imgName '_x' num2str(x) '_y' num2str(y) '.mat'];

%%-- Option 2: direct from image file

% maskNucFile='example_data/nuclei_mask_1.png';
% maskEpiFile='example_data/epithelium_mask_1.png';
% lympDataFile='example_data/lymp_data_1.mat';
% I=imread('example_data/image_1.png');

%% Showing original image

%imshow(I); % Shows the tile

%% Reading nuclei mask

MN=imread(maskNucFile);
if size(MN,3)>1
    MN(:,:,1:2)=[];
end
MN=imbinarize(MN,graythresh(MN));

%imshow(imoverlay(I,MN,[0 1 0])); % Overlays nuclei (no alpha)
%imshow(I); alphamask(logical(MN),[0 0 1],.1); % Overlays nuclei (with alpha)

%% Reading epithelium mask

ME=imread(maskEpiFile);
if size(ME,3)>1
    ME=rgb2gray(ME);
end
%ME=imresize(ME,4);

%imshow(I); alphamask(logical(ME),[0 0 1],.3); % Overlays epithelium (with alpha)

%% Reading lymphocyte-related data
data=load(lympDataFile);
classes=getEpiStrLympClasses(data.isLymphocyte,~data.epiNuc);
numCent=length(data.nucCentroids);

%drawNucContoursByClass(MN,I,data.nucCentroids,ones(numCent,1)*2,5); % Plots nuclei borders using a single color
%drawNucleiCentroidsByClass(I,data.nucCentroids,ones(numCent,1)*3,12); % Plots nuclei as dots using a single color

%drawNucleiCentroidsByClass(I,data.nucCentroids,classes,10); % Plots nuclei as dots with diff. color per class.
%drawNucContoursByClass(MN,I,data.nucCentroids,classes,5); % Plots nuclei borders with diff. color per class.

%% Cell cluster graphs

coords= {   data.nucCentroids(~data.isLymphocyte & ~data.epiNuc,:),...
            data.nucCentroids(data.isLymphocyte & data.epiNuc,:),...
            data.nucCentroids(data.isLymphocyte & ~data.epiNuc,:),...            
            data.nucCentroids(~data.isLymphocyte & data.epiNuc,:),...
        };

numGroups=length(coords);
alpha=.43;
r=.185;
M=cell(numGroups,1);
for i=1:numGroups
    [~,~,~,~,groupMatrix] = construct_nodesCluster(struct('centroid_r',coords{i}(:,2)','centroid_c',coords{i}(:,1)'), alpha, r);
    M{i}=groupMatrix;
end

%drawGraph(I,coords,M);
%drawGraph(I,coords,M,5,20,.8,.8,false);

%args.printGroupId=true;
args.fillAlpha=.3;
args.lineWidth=3;
args.lineAlpha=.7;
args.useBlackLine=true;
%args.plotPoints=false;
%args.printCenter=true;
%args.plotCenterGrouping=true;

plotClusters(I,coords,M,args);

