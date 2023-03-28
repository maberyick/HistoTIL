I=imread('data/case713_1001_501.png'); % original image
M=imread('data/case713_1001_501_mask.png'); % mask (nuclei segmentation)

% If mask is not available, use:
%M=getWatershedMask(I);

lympModel=load('data/model.mat'); % trained models

% centroids contains the position (x,y) of all detected nuclei
% isLympocyte indicates wheter or not each detected nucleus is a lymphocyte. (1 - lymphocyte, 0 - non-lymphocyte)
[centroids,isLymphocyte]=runLympClassificationModel(I,M,lympModel.model);

% Optional: draws the image showing the class of each centroid (lymphocyte or non-lymphocyte)
drawNucleiCentroidsByClass(I,centroids,isLymphocyte);
