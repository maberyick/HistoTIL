function [nucleiCentroids,prediction]=runLympClassificationModel(Image,Mask,model )
%RUNLYMPCLASSIFICATIONMODEL runs a classification model that separates
%lymphocytes (1) from non-lymphocytes (2) for an input image
% Image - Input RGB image
% Maks - Mask image(binary or grayscale)
% model - The SVM separation model

% Extraction nuclei features
[nucleiCentroids,nucleiFeatures,~]=getNucleiFeatures(Image,Mask);

% Using the SVM model
numNuc=size(nucleiFeatures,1);

% using libsvm
%[prediction,~,~] = svmpredict(ones(numNuc,1), nucleiFeatures, model,'-q');

%using matlab
prediction = predict(model,nucleiFeatures);


end

