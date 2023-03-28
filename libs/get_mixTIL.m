function [nucleiCentroids,isLymphocyte,nucFeatures,denFeat,spaFeat,ctxFeat,spaFeat_v2,perinuclear,nucleiCentroids_stro,isLymphocyte_stro,nucFeatures_stro,denFeat_stro,spaFeat_stro,ctxFeat_stro,spaFeat_v2_stro,perinuclear_stro,nucleiCentroids_epi,isLymphocyte_epi,nucFeatures_epi,denFeat_epi,spaFeat_epi,ctxFeat_epi,spaFeat_v2_epi,perinuclear_epi,nucleiCentroids_bund,isLymphocyte_bund,nucFeatures_bund,denFeat_bund,spaFeat_bund,ctxFeat_bund,spaFeat_v2_bund,perinuclear_bund] = get_mixTIL(I,ES,M,maskFile,maskEFile,maskSFile,maskBFile,ESmaskFile,SSmaskFile,EBmaskFile,lympModel)
%GETPERINUCLEARFEATURES
% Clean image by watershed
% Check if the mask is already processed
if exist(maskFile,'file')~=2
    [M] = get_cleanMask(I,M,maskFile);
    M = imbinarize(M);
else
    % load mask
    M = imread(maskFile);
    M = imbinarize(M);
end
%%Adapt the ES mask to the nuclei
if [exist(maskEFile,'file') && exist(maskSFile,'file') && exist(maskBFile,'file')] == 0
    [M_stro,M_epi,M_bund] = get_epistroma(I,M,ES,maskEFile,maskSFile,maskBFile,ESmaskFile,SSmaskFile,EBmaskFile);
else
    % load mask
    M_epi = imread(maskEFile);
    M_stro = imread(maskSFile);
    M_bund = imread(maskBFile);
end
% Extract based on all tissue cells
[nucleiCentroids,nucFeatures] = getNucLocalFeatures(I,M);
if sum(sum(M)) == 0
    spaFeat = zeros(1,85);
    nucFeatures = zeros(1,100);
    denFeat = zeros(1,19);
    ctxFeat = zeros(1,87);
    spaFeat_v2 = zeros(1,1400);
else
    isLymphocyte = (predict(lympModel,nucFeatures(:,1:7)))==1;
    lympCentroids=nucleiCentroids(isLymphocyte==1,:);
    nonLympCentroids=nucleiCentroids(isLymphocyte~=1,:);
    nucAreas=nucFeatures(:,1);
    lympAreas=nucAreas(isLymphocyte==1);
    %% SpaTIL (v1) - 85
    [spaFeat]=getSpaTILFeatures(lympCentroids,nonLympCentroids);
    %% denTIL (v1) - 19
    [denFeat]=getDenTILFeatures(I,lympCentroids,nonLympCentroids,lympAreas);
    %% contextTIL (v1) - 87
    [ctxFeat]=getContextTILFeatures(nucleiCentroids,nucAreas,nucFeatures(:,2),nucFeatures(:,8),nucFeatures(:,11),nucFeatures(:,5),nucFeatures(:,13),nucFeatures(:,14));
    %% nucFeat (v1) - 100
    % [nucleiCentroids,nucFeatures] = getNucLocalFeatures(I,M);
    %% spaTIL (v2) - 350
    [spaFeat_v2]=getSpaTILFeaturesv2(lympCentroids,nonLympCentroids);
    %% peri-nuclear Features
    [perinuclear]=getPeriNuclearFeatures(I,M);
end
% Extract based on the Epistroma cells
[nucleiCentroids_epi,nucFeatures_epi] = getNucLocalFeatures(I,M_epi);
if sum(sum(M_epi)) == 0
    spaFeat_epi = zeros(1,85);
    nucFeatures_epi = zeros(1,100);
    denFeat_epi = zeros(1,19);
    ctxFeat_epi = zeros(1,87);
    isLymphocyte_epi = 0;
else
    try
        isLymphocyte_epi = (predict(lympModel,nucFeatures_epi(:,1:7)))==1;
    catch
        isLymphocyte_epi = 0;
    end
    if isLymphocyte_epi == 0
        nucFeatures_epi = zeros(1,100);
        spaFeat_epi = zeros(1,85);
        denFeat_epi = zeros(1,19);
        ctxFeat_epi = zeros(1,87);
    else
        lympCentroids_epi=nucleiCentroids_epi(isLymphocyte_epi==1,:);
        nonLympCentroids_epi=nucleiCentroids_epi(isLymphocyte_epi~=1,:);
        nucAreas_epi=nucFeatures_epi(:,1);
        lympAreas_epi=nucAreas_epi(isLymphocyte_epi==1);
        [spaFeat_epi]=getSpaTILFeatures(lympCentroids_epi,nonLympCentroids_epi);
        [denFeat_epi]=getDenTILFeatures(I,lympCentroids_epi,nonLympCentroids_epi,lympAreas_epi);
        [ctxFeat_epi]=getContextTILFeatures(nucleiCentroids_epi,nucAreas_epi,nucFeatures_epi(:,2),nucFeatures_epi(:,8),nucFeatures_epi(:,11),nucFeatures_epi(:,5),nucFeatures_epi(:,13),nucFeatures_epi(:,14));
        [spaFeat_v2_epi]=getSpaTILFeaturesv2(nucleiCentroids_epi,nonLympCentroids_epi);
        [perinuclear_epi]=getPeriNuclearFeatures(I,M_epi);
    end
end
% Extract based on the Stroma cells
[nucleiCentroids_stro,nucFeatures_stro] = getNucLocalFeatures(I,M_stro);
if sum(sum(M_stro)) == 0
    spaFeat_stro = zeros(1,85);
    nucFeatures_stro = zeros(1,100);
    denFeat_stro = zeros(1,19);
    ctxFeat_stro = zeros(1,87);
    isLymphocyte_stro = 0;
else
    try
        isLymphocyte_stro = (predict(lympModel,nucFeatures_stro(:,1:7)))==1;
    catch
        isLymphocyte_stro = 0;
    end
    if isLymphocyte_stro == 0
        spaFeat_stro = zeros(1,85);
        nucFeatures_stro = zeros(1,100);
        denFeat_stro = zeros(1,19);
        ctxFeat_stro = zeros(1,87);
    else
        lympCentroids_stro=nucleiCentroids_stro(isLymphocyte_stro==1,:);
        nonLympCentroids_stro=nucleiCentroids_stro(isLymphocyte_stro~=1,:);
        nucAreas_stro=nucFeatures_stro(:,1);
        lympAreas_stro=nucAreas_stro(isLymphocyte_stro==1);
        [spaFeat_stro]=getSpaTILFeatures(lympCentroids_stro,nonLympCentroids_stro);
        [denFeat_stro]=getDenTILFeatures(I,lympCentroids_stro,nonLympCentroids_stro,lympAreas_stro);
        [ctxFeat_stro]=getContextTILFeatures(nucleiCentroids_stro,nucAreas_stro,nucFeatures_stro(:,2),nucFeatures_stro(:,8),nucFeatures_stro(:,11),nucFeatures_stro(:,5),nucFeatures_stro(:,13),nucFeatures_stro(:,14));
        [spaFeat_v2_stro]=getSpaTILFeaturesv2(nucleiCentroids_stro,nonLympCentroids_stro);
        [perinuclear_stro]=getPeriNuclearFeatures(I,M_stro);
    end
end
% Extract based on the Stroma cells
[nucleiCentroids_bund,nucFeatures_bund] = getNucLocalFeatures(I,M_bund);
if sum(sum(M_bund)) == 0
    spaFeat_bund = zeros(1,85);
    nucFeatures_bund = zeros(1,100);
    denFeat_bund = zeros(1,19);
    ctxFeat_bund = zeros(1,87);
    isLymphocyte_bund = 0;
else
    try
        isLymphocyte_bund = (predict(lympModel,nucFeatures_bund(:,1:7)))==1;
    catch
        isLymphocyte_bund = 0;
    end
    if isLymphocyte_bund == 0
        spaFeat_bund = zeros(1,85);
        nucFeatures_bund = zeros(1,100);
        denFeat_bund = zeros(1,19);
        ctxFeat_bund = zeros(1,87);
    else
        lympCentroids_bund=nucleiCentroids_bund(isLymphocyte_bund==1,:);
        nonLympCentroids_bund=nucleiCentroids_bund(isLymphocyte_bund~=1,:);
        nucAreas_bund=nucFeatures_bund(:,1);
        lympAreas_bund=nucAreas_bund(isLymphocyte_bund==1);
        [spaFeat_bund]=getSpaTILFeatures(lympCentroids_bund,nonLympCentroids_bund);
        [denFeat_bund]=getDenTILFeatures(I,lympCentroids_bund,nonLympCentroids_bund,lympAreas_bund);
        [ctxFeat_bund]=getContextTILFeatures(nucleiCentroids_bund,nucAreas_bund,nucFeatures_bund(:,2),nucFeatures_bund(:,8),nucFeatures_bund(:,11),nucFeatures_bund(:,5),nucFeatures_bund(:,13),nucFeatures_bund(:,14));
        [spaFeat_v2_bund]=getSpaTILFeaturesv2(nucleiCentroids_bund,nonLympCentroids_bund);
        [perinuclear_bund]=getPeriNuclearFeatures(I,M_bund);
    end
end
end