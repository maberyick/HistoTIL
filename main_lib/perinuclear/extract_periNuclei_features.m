function [feats_periNuclei,description_periNuclei] =  extract_periNuclei_features(img, img2, cellRadius, Area, binfo)
%cytoplasm area (periNuclei)
%the ratio of cytoplasm area (periNuclei) to nuclei area
%haralick features in periNuclei region
%entropy of the perinuclei region



%Notice: the haralick input is the masked rgb images, so it's necessary to
%generate the mask of peri-nuclei region. Needs to treat the area and
%texture feature separately, since the texture will calculate based on the
%co-occurence matrix of entire matrix. The area should calculated based per
%object

periAreas = pi*cellRadius^2-Area;
areaRatioCyto2Nuclei = periAreas./Area;


%calculating haralick feature
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
grays = info.grays;
win = info.win;
dist = info.dist;


himg = uint16(rescale_range(gimg,0,grays-1));
%% Adapted code for the new haralick
centroids = cat(1, binfo.Centroid);
harFeat=[];
bbx_rad = 5;
for i=1:length(centroids)
    nucleus=binfo(i);
    bbox = nucleus.BoundingBox;
    bbox = [round(bbox(1)) round(bbox(2)) (bbox(3) - 1) (bbox(4) - 1)];
    try
        roi = img2(bbox(2)-bbx_rad : bbox(2) + bbox(4)+bbx_rad, bbox(1)-bbx_rad : bbox(1) + bbox(3)+bbx_rad, :);
    catch merror
        if ~isempty(merror)
            roi = img2(bbox(2) : bbox(2) + bbox(4), bbox(1) : bbox(1) + bbox(3), :);
        end
    end
    merror = [];
    %imshow(roi)
    gray=rgb2gray(roi);
    % Haralick features
    glcm = graycomatrix(gray);    
    harFeat=[harFeat;haralickTextureFeatures(glcm,1:14)'];
end
%%

Texturefeats = [mean(harFeat) std(harFeat)];

J =  entropy(img);

feats_periNuclei = [mean(periAreas), median(periAreas), var(periAreas), std(periAreas),...
    mean(areaRatioCyto2Nuclei), median(areaRatioCyto2Nuclei), var(areaRatioCyto2Nuclei), std(areaRatioCyto2Nuclei),...
    Texturefeats, J];
description_periNuclei = {'mean value of periNuclei area', 'median of periNuclei area', 'variance of periNuclei area', 'standard deviation of periNuclei area', ...
    'mean value of periNuclei area to nuclei area ratio', 'mean value of periNuclei area to nuclei area ratio', 'mean value of periNuclei area to nuclei area ratio', 'mean value of periNuclei area to nuclei area ratio',...
    'Haralick:mean intensity angular_second_moment_energy','Haralick:standard deviation intensity angular_second_moment_energy',...
    'Haralick:mean intensity contrast_energy','Haralick:standard deviation intensity contrast_energy',...
    'Haralick:mean intensity correlation','Haralick:standard deviation intensity correlation',...
    'Haralick:mean intensity contrast_var','Haralick:standard deviation intensity contrast_var',...
    'Haralick:mean intensity inverse_difference_moment_homogeneity','Haralick:standard deviation intensity inverse_difference_moment_homogeneity',...
    'Haralick:mean intensity sum_average','Haralick:standard deviation intensity sum_average',...
    'Haralick:mean intensity sum_variance_aprox','Haralick:standard deviation intensity sum_variance_aprox',...
    'Haralick:mean intensity sum_entropy','Haralick:standard deviation intensity sum_entropy',...
    'Haralick:mean intensity entropy_cut_out_zeros','Haralick:standard deviation intensity entropy_cut_out_zeros',...
    'Haralick:mean intensity difference_variane,','Haralick:standard deviation intensity difference_variane,',...
    'Haralick:mean intensity difference_entropy','Haralick:standard deviation intensity difference_entropy',...
    'Haralick:mean intensity information_measure_of_correlation_I','Haralick:standard deviation intensity information_measure_of_correlation_I',...
    'Haralick:mean intensity information_measure_of_correlation_II','Haralick:standard deviation intensity information_measure_of_correlation_II',...
    'Haralick:mean intensity maximal_correlation_coefficient','Haralick:standard deviation intensity maximal_correlation_coefficient',...
    'entropy of periNuclei pixel intensity'};
end
