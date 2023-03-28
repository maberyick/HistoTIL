clc; clear all;
path = '/media/maberyick/MABERYICK1/dataset/ETMA_Chemo_Set/test_sample/A004_J10.mat';
load(path)
load([path(1:end-4),'_mask.mat'])
HEtile = tileStruct(1).data;
nucleiMaskTile = maskStruct(1).nucmsk;
nucleiMaskTile = logical(nucleiMaskTile);

HEtile = HEtile(501:1500,501:1500,:);
nucleiMaskTile = nucleiMaskTile(501:1500,501:1500,:);
        
%overlaid = imoverlay(HEtile,nucleiMaskTile,[0 1 0]);
%imshow(overlaid)
%%
addpath(genpath('/home/maberyick/Dropbox/main_ResDev/tool_featureExtraction/pj_periNuclear_V2/FeatureExtraction/'))

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
img = HEtile.*maskPeriNuclei;
img2 = HEtile.*uint8(~nucleiMaskTile);
Area =  [binfo.ConvexArea];
%[feats_periNuclei,description_periNuclei] = extract_periNuclei_features(img1, img2, cellRadius, Area, binfo);

%
%% texture
% alternatively 13 haralick
info.dist = 1;
info.win = 1;
info.grays = 256;

n = 0;

if ndims(img) < 3
    gimg = img;
elseif ndims(img) == 3
    gimg = rgb2gray(img); % assume img is rgb
else
    fprintf('Unidentified image format')
end

mask(:,:,1) = (gimg ~= max(max(gimg)));
%for grays = [64,128,256]
%    for win = [1,2,3]
%        n = n + 1
grays = info.grays;
win = info.win;
dist = info.dist;

himg = uint16(rescale_range(gimg,0,grays-1));
f = haralick_img(himg,mask,grays,win,dist,1);
Hf = f.img3;
%%
imshow(HEtile)
figure
%hold on
mask2_har = Hf(:,:,2);
imshow(mask2_har)

test = mask2_har;
%test(test < 0.01)=0;

thresh = multithresh(test ,8);
seg_I = imquantize(test ,thresh);
RGB = label2rgb(seg_I);
figure;
imshow(test)
figure;
imshow(RGB)
%colormap gray
%RGB(:,:,3) = 255;  % Remove all blue components from every pixel.
%imshow(RGB)


%colormap('summer')
