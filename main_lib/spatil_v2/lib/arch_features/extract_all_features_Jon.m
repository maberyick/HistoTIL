function [feats,description] = extract_all_features_Jon(mask,img,descriptor, varagin)
%% Extracts features from nuclei/glandular bounds
% bounds is a struct containing the centroid and boundary info
% img is a color or grayscale image used for computed for haralick features
% img can be omitted if not using haralick featuers

% if i really want to fix those names
% if strcmp(type, 'nuclei')
%      load NewImgDescription.mat
%     temp = regexp(o.Description.ImageFeatures,'CGT');
%     temp2 = cellfun(@isempty,temp,'UniformOutput',false);
%
%     %description{3} = o.Description.ImageFeatures([temp2{:}] == 0);
% else
    load NewImgDescription2.mat
% end

%% debug
    %bounds=epbounds
%%


% Dimred2 = @(x) Dimred(x, 'GE', param);
% EigV = cellfun(Dimred2,data,'UniformOutput',false);

% if nargin < 3
%     mask(:,:,1) = (img ~= max(max(img)));
% end
if nargin < 3
    descriptor = '';
end
bounds = DLMask2bounds_function(img, mask);
% [stbounds] = DLMask2bounds_function(img, stromamask);
% mask=bounds2mask(bounds, img);

if ischar(img)
    img=imread(img);
end

if size(img,3)==4
    img=img(:,:,1:3);
end
if ischar(mask)
    mask=imread(mask);
end
fprintf('\nExtracting Graph Features...')
% try
feats{1} = extract_graph_feats(bounds);
% catch
%    feats{1}=NaN(1,51);
% end
description{1} = [o.Description.ImageFeatures(1:51)];

% fprintf('\nExtracting Morph Features...')
% try
feats{2} = extract_morph_feats(bounds);
% catch
%     feats{2}=NaN(1,100);
% end
temp = regexp(o.Description.ImageFeatures,'Shape');
temp2 = cellfun(@isempty,temp,'UniformOutput',false);
description{2} = [o.Description.ImageFeatures([temp2{:}] == 0)];

% fprintf('\nExtracting CGT Features...')
%% CGT feats
%These Cell Graph Tensors are George's orientation graphs
%It measures the entropy of the nuclei/gland orientations

% try
feats{3} = extract_CGT_feats(bounds);
% catch
%     feats{3}=NaN(1,39);
% end
temp = regexp(o.Description.ImageFeatures,'CORE');
temp2 = cellfun(@isempty,temp,'UniformOutput',false);
description{3} = [o.Description.ImageFeatures([temp2{:}] == 0)];

%  fprintf('\nExtracting Cluster Graph Features...')
% try
%The GSG features are the features having to do with the neighborhoods of 
%nuclei or glands, and how they are spaced in relationship to each other
%% Cluster graph feats
    [tempclusterfeats,~]=extract_cluster_graph_feats(bounds);
    tempclusterfeats=cellvec2array(tempclusterfeats);
feats{4} = tempclusterfeats;
% catch
%     feats{4}=NaN(1,26);
% end
temp = regexp(o.Description.ImageFeatures,'CCG');
temp2 = cellfun(@isempty,temp,'UniformOutput',false);
description{4} = [o.Description.ImageFeatures([temp2{:}] == 0)];

if nargin > 1
%     fprintf('\nExtracting Texture Features...')
%     try
%     feats{5} = extract_texture_feats(img, mask);
% %     catch
% %        feats{5}=NaN(1,26);
% %     end
%     temp = regexp(o.Description.ImageFeatures,'Haralick');
%     temp2 = cellfun(@isempty,temp,'UniformOutput',false);
%     description{5} = [o.Description.ImageFeatures([temp2{:}] == 0)];
    
    
    description=[description{:}];
    for i=1:length(description)
        description{i}=[description{i}, descriptor];
    end
end

%feat_description


function [graphfeats] = extract_graph_feats(bounds)
gb_r = [bounds.centroid_r];
gb_c = [bounds.centroid_c];

if length(bounds) >=3 %Delaunay triangle requires at least 3 points
    [graphfeats] = get_graph_features([gb_r]',[gb_c]');
else
    [graphfeats]=NaN(1,51);
end


function [morphfeats]=extract_morph_feats(bounds)
%% Morph
gb_r = {bounds.r};
gb_c = {bounds.c};

badglands = [];
for j = 1:length(gb_r)
    try
        [feat] = morph_features([gb_r{j}]',[gb_c{j}]');
        feats(j,:) = feat;
        
    catch ME
        badglands = [badglands j]; %track bad glands
    end
end

feats(badglands,:) = []; %remove bad glands

morphfeats = [mean(feats) std(feats) median(feats) min(feats)./max(feats)];



function [CGTfeats] = extract_CGT_feats(bounds)
%% CGT
% addpath(genpath(['Z:\2012_10_GlandOrientation']))

%[CGTfeats,c_matrix,info] = get_CGT_features(bounds);
a=0.5;
r=0.2;
%[CGTfeats, c_matrix, info, feats, network, edges] = get_CGT_features_networks_weight(bounds,a,r);
if length(bounds) >=3
    [CGTfeats, c_matrix, info, feats, network, edges] = extract_CGT_features(bounds,a,r);
else
    CGTfeats=NaN(1,39);
end

function [CCGfeats,feature_list] = extract_cluster_graph_feats(bounds)
%CCG

%addpath(genpath(['Z:\2012_10_GlandOrientation']))

info.alpha = 0.5;
info.radius = 0.2;
% build graph
% alpha = info.alpha;
% r = info.radius;

alpha = .45; %These are alpha and r values I picked based on how they looked on a few images.  Very arbitrary, but
%small changes make a huge difference.  Will need to look into this
%further.
r = .21;
%[VX,VY,x,y,edges] = construct_ccgs(bounds,alpha, r);
[VX,VY,x,y,edges] = construct_ccgs_optimized(bounds,alpha, r);

%[CCGfeats,feature_list] = cluster_graph_features_networkboost(bounds, edges);
if length(bounds) >=3
    [CCGfeats,feature_list] = cluster_graph_features_optimized(bounds, edges);
else
    CCGfeats=NaN(1,26);
    feature_list=[];
end

function [Texturefeats] = extract_texture_feats(img, mask)

%% texture
%addpath(genpath(['Z:\2012_10_GlandOrientation']))

% alternatively 13 haralick
%addpath(genpath('Z:\Datasets\SatishData\haralick'))

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

if ~exist('mask', 'var');
    mask(:,:,1) = (gimg ~= max(max(gimg)));
end
%for grays = [64,128,256]
%    for win = [1,2,3]
%        n = n + 1
grays = info.grays;
win = info.win;
dist = info.dist;

himg = uint16(rescale_range(gimg,0,grays-1));
%himg = rescale_range(gimg,0,grays-1);
%f = haralick_img(himg,mask,grays,win,dist,1);
f = haralick_img(himg,mask,grays,win,dist,1);
Hf = f.img3;
%        HaralickFeat{n} = single(Hf);
%    end
%end

for i = 1:size(Hf,3)
    feat = Hf(:,:,i);
    img_mask = mask(:,:,1);
    roi = feat(img_mask ==1);
    
    Texturefeats(n+1) = mean(roi(:));
    Texturefeats(n+2) = std(roi(:));
    %Texturefeats(n+3) = mode(roi(:));
    n = n + 2;
end


% count = 1;
% modifier = [{'mean '} {'standard deviation '}]
% for j = 1:numel(feat.names)
% for i = 1:numel(modifier)
%     HaralickTextureFeatureDescription{count} = [modifier{i} 'intensity ' feat.names{j}];
%     count = count + 1;
% end
% end
