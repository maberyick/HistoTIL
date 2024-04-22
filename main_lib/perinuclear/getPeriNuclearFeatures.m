function [features] = getPeriNuclearFeatures(HEtile,nucleiMaskTile)
%GETPERINUCLEARFEATURES Gets a set of features from the peri-nuclear (membrane of the cell)

%         HEtile = HEtile(51:100,51:100,:);
%         nucleiMaskTile = nucleiMaskTile(51:100,51:100,:);

nucleiMaskTile = logical(nucleiMaskTile);
%overlaid = imoverlay(HEtile,nucleiMaskTile,[0 1 0]);
%imshow(overlaid)
%close all

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
%img = HEtile.*maskPeriNuclei;
%Area = [binfo.ConvexArea];

%mask2 = logical(maskPeriNuclei(:,:,1));

[features,~] = extract_periNuclei_features(HEtile.*maskPeriNuclei, HEtile.*uint8(~nucleiMaskTile), cellRadius, [binfo.ConvexArea], binfo);
%featPerPatient = [featPerPatient;feats_periNuclei];
end