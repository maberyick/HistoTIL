clear,clc, close all
if isunix
    imageTilesPath = '/mnt/projects/CSE_BME_AXM788/data/CCF/Lung/subset_ETMA/mat_tiles/';
    maskTilesPath = '/mnt/projects/CSE_BME_AXM788/data/CCF/Lung/subset_ETMA/mat_nucmasks/';
    savingPath = '/scratch/users/xxw345/';
    addpath(genpath('/home/xxw345/matlabWorkSpace/dependencies'))
elseif ispc
    imageTilesPath = 'Z:\data\CCF\Lung\subset_ETMA\mat_tiles\';
    maskTilesPath = 'Z:\data\CCF\Lung\subset_ETMA\mat_nucmasks\';
    addpath(genpath('C:\Users\xxw345\codes\FeatureExtraction\'))
end
%%
clc; clear all;
savingPath = '/media/maberyick/MABERYICK1/dataset/testing_folder/tile/saved/';
imageTilesPath = '/media/maberyick/MABERYICK1/dataset/testing_folder/tile/';
maskTilesPath = '/media/maberyick/MABERYICK1/dataset/testing_folder/mask/';
%%
imageTiles = dir([imageTilesPath,'*.mat']);
maskTiles = dir([maskTilesPath,'*.mat']);

for i = 1:length(imageTiles)
    imageTiles(i).name
    load([imageTilesPath,imageTiles(i).name])
    load([maskTilesPath,imageTiles(i).name(1:end-4) '_mask.mat'])
    
    if exist([savingPath,imageTiles(i).name],'file') == 2
        continue
    end
%     save([savingPath,imageTiles(i).name],'imageTiles(i).name')
    featPerPatient = [];
    for j = 1:length(tileStruct)
        HEtile = tileStruct(j).data;
        nucleiMaskTile = maskStruct(j).nucmsk;
        
%         HEtile = HEtile(51:100,51:100,:);
%         nucleiMaskTile = nucleiMaskTile(51:100,51:100,:);
        
        nucleiMaskTile = logical(nucleiMaskTile);                
        
        overlaid = imoverlay(HEtile,nucleiMaskTile,[0 1 0]);
        imshow(overlaid)
        close all
        
        maskperim = bwperim(nucleiMaskTile);
        binfo = regionprops(maskperim,'ConvexHull', 'Centroid', 'ConvexArea', 'BoundingBox');
        centers = cat(1,binfo.Centroid);
        cellRadius = 15;
        if size(centers,1) == 0
            mask = zeros(size(nucleiMaskTile));
        else
            mask = createCirclesMask(size(nucleiMaskTile), centers, cellRadius*ones(size(centers,1),1));
        end
        maskPeriNuclei = mask-nucleiMaskTile;
        maskPeriNuclei = uint8(repmat(maskPeriNuclei,1,1,3));
        % Test it
        img = HEtile.*maskPeriNuclei;
        Area = [binfo.ConvexArea];
        
        mask2 = logical(maskPeriNuclei(:,:,1));

        [feats_periNuclei,description_periNuclei] = extract_periNuclei_features(HEtile.*maskPeriNuclei, HEtile.*uint8(~nucleiMaskTile), cellRadius, [binfo.ConvexArea], binfo);
        featPerPatient = [featPerPatient;feats_periNuclei];
         nuclei = {binfo.ConvexHull};
         bounds = veta2bounds(nuclei);
%         nucleiMaskTile = uint8(repmat(nucleiMaskTile,1,1,3));
%         [feats,description] = extract_all_features(bounds,HEtile.*nucleiMaskTile);
         [feats,description] = extract_all_features(bounds,nucleiMaskTile,HEtile);
         featPerPatient = [featPerPatient;feats];
    end
         save([savingPath,imageTiles(i).name],'featPerPatient','description_periNuclei')
    clear featPerPatient
end