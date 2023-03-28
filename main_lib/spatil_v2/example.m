addpath(genpath('/media/maberyick/Data4/CCIPD_data/NSCLC_Chemo/mat_libraries/SpaTIL/'))
clear;clc;

I=imread('/media/maberyick/Data4/CCIPD_data/NSCLC_Chemo/mat_libraries/SpaTIL/example_data/image_1.png');
M=imread('/media/maberyick/Data4/CCIPD_data/NSCLC_Chemo/mat_libraries/SpaTIL/example_data/nuclei_mask_1.png');

[nucleiCentroids,nucFeatures] = getNucLocalFeatures(I,M);

model=load('/media/maberyick/Data4/CCIPD_data/NSCLC_Chemo/mat_libraries/SpaTIL/example_data/lymp_svm_matlab_40x.mat');
isLymphocyte = (predict(model.model,nucFeatures(:,1:7)))==1;

coords={nucleiCentroids(isLymphocyte,:),nucleiCentroids(~isLymphocyte,:),};

figure;imshow(I);
hold on;
%[feat,featNames]=getSpaTILFeatures_v2(coords,5,.44);
[feat,featNames]=getSpaTILFeatures_v2(coords);
hold off;
