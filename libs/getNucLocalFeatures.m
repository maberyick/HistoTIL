function [ centroids,features] = getNucLocalFeatures( image, mask )
%GETNUCLOCALFEATURES Summary of this function goes here
%   Detailed explanation goes here
orderZernikePol=7;
image=normalizeStaining(image);
grayImg=rgb2gray(image);
regionProperties = regionprops(mask,grayImg,'Centroid','Area',...
    'BoundingBox','Eccentricity','EquivDiameter','Image',...
    'MajorAxisLength','MaxIntensity','MeanIntensity','MinIntensity',...
    'MinorAxisLength','Orientation','PixelValues');
centroids = cat(1, regionProperties.Centroid);
nucleiNum=size(regionProperties,1);
medRed=zeros(nucleiNum,1);
medGreen=zeros(nucleiNum,1);
medBlue=zeros(nucleiNum,1);
entropyRed=zeros(nucleiNum,1);
medIntensity=zeros(nucleiNum,1);
entropyIntensity=zeros(nucleiNum,1);
edgeMedIntensity=zeros(nucleiNum,1);
harFeat=zeros(nucleiNum,14);
zernFeat=zeros(nucleiNum,36*2);
parfor i=1:nucleiNum
    nucleus=regionProperties(i);
    bbox = nucleus.BoundingBox;
    bbox = [round(bbox(1)) round(bbox(2)) (bbox(3) - 1) (bbox(4) - 1)];
    roi = parimgbox_mix(image, bbox);
    per = bwperim(nucleus.Image);
    gray=rgb2gray(roi);
    R=roi(:,:,1);
    G=roi(:,:,2);
    B=roi(:,:,3); 
    R=R(nucleus.Image == 1);
    G=G(nucleus.Image == 1);
    B=B(nucleus.Image == 1);
    grayPix=gray(nucleus.Image == 1);
    perPix=gray(per==1);  
    % Intensity features:
    medRed(i) = median(double(R));
    medGreen(i) = median(double(G));
    medBlue(i) = median(double(B));
    medIntensity(i)=median(double(grayPix));
    edgeMedIntensity(i)=median(double(perPix));
    % Texture features:
    % Entropies
    entropyRed(i)=getNucEntropy(R);
    entropyIntensity(i)=getNucEntropy(grayPix);
    % Haralick features
    glcm = graycomatrix(gray);
    harFeat(i,:)=haralickTextureFeatures(glcm,1:14)';
    % Zernike moments
    [w,h]=size(nucleus.Image);
    sqRoi=getSquareRoi(grayImg,nucleus.Centroid,w,h);
    zernFeat(i,:)=getZernikeMomentsImg(sqRoi,orderZernikePol);
end
ratioMedRB=medRed./medBlue;
ratioMedRG=medRed./medGreen;
ratioAxes=[regionProperties.MajorAxisLength]./[regionProperties.MinorAxisLength];
features = horzcat([regionProperties.Area]',[regionProperties.Eccentricity]', ...
    ratioAxes',medRed,entropyRed,double([regionProperties.MinIntensity]'),...
    double([regionProperties.MaxIntensity]'),[regionProperties.EquivDiameter]',...
    [regionProperties.Orientation]',entropyIntensity,medIntensity,edgeMedIntensity,...
    ratioMedRB,ratioMedRG,harFeat,zernFeat...
    );
end

