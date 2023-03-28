% Get the spatil (350 features) from the patches
% Linux
%addpath(genpath('/media/maberyick/Data4/CCIPD_data/NSCLC_Chemo/mat_libraries/'))
%addpath(genpath('/media/maberyick/Data4/CCIPD_data/NSCLC_Chemo/'))
% Windows
addpath(genpath('E:\CCIPD_Projects\PhenoTIL_Chemo\ML_approach\mat_libraries\'))
addpath(genpath('E:\CCIPD_Projects\PhenoTIL_Chemo\ML_approach\'))
% image path
% Linux
%img_path = '/media/maberyick/Data4/CCIPD_data/NSCLC_Chemo/ccf_cohort/patch_images/';
%mask_lymp_path = '/media/maberyick/Data4/CCIPD_data/NSCLC_Chemo/ccf_cohort/masks/lymp/';
%mask_nonlymp_path = '/media/maberyick/Data4/CCIPD_data/NSCLC_Chemo/ccf_cohort/masks/nonlymp/';
%spatil_feature_path = '/media/maberyick/Data4/CCIPD_data/NSCLC_Chemo/ccf_cohort/features/spatil2/';
% Windows
img_path = 'E:\CCIPD_Projects\PhenoTIL_Chemo\ML_approach\bms_cohort\patch_images\';
mask_lymp_path = 'E:\CCIPD_Projects\PhenoTIL_Chemo\ML_approach\bms_cohort\masks\lymp\';
mask_nonlymp_path = 'E:\CCIPD_Projects\PhenoTIL_Chemo\ML_approach\bms_cohort\masks\nonlymp\';
spatil_feature_path = 'E:\CCIPD_Projects\PhenoTIL_Chemo\ML_approach\bms_cohort\features\spatil2\';
% Get the dir of files
img_dir = dir([img_path '*.png']);
% Linux
%lympModel = load('./mat_libraries/lymp_svm_matlab_wsi.mat');
% Windows
lympModel = load('E:\CCIPD_Projects\PhenoTIL_Chemo\ML_approach\mat_libraries\lymp_svm_matlab_wsi.mat');
% loop through the files
% loop through the values

% h = waitbar(0,'Computing spaTIL...');
pb = CmdLineProgressBar('Computing spaTIL...');

for x=1:length(img_dir)
    %perct = x/length(img_dir);
    patch_name = img_dir(x).name;
    %fprintf('\n')
    %disp(patch_name)
    % Check if files exists
    if isfile([spatil_feature_path patch_name(1:end-4) '.mat'])
        %disp('feature exists')
        continue
    else
        %load the image
        patch_image = imread([img_path patch_name]);
        try
            % load the lymp mask
            lymp_mask = imread([mask_lymp_path patch_name]);
            lymp_mask = logical(lymp_mask);
            % load the non lymp mask
            nonlymp_mask = imread([mask_nonlymp_path patch_name]);
            nonlymp_mask = logical(nonlymp_mask);
        catch
            %disp('mask missing; calculating...')
            % calculate the mask
            patch_mask = logical(getWatershedMask(patch_image,true,6,10));
            [Centroids,nucFeatures] = getNucLocalSimpleFeatures(patch_image,patch_mask);
            isLymphocyte = (predict(lympModel.model,nucFeatures(:,1:7)))==1;
            lympCentroids=Centroids(isLymphocyte==1,:);
            nonLympCentroids=Centroids(isLymphocyte~=1,:);
            lymp_mask = bwselect(patch_mask,lympCentroids(:,1),lympCentroids(:,2));
            nonlymp_mask = bwselect(patch_mask,nonLympCentroids(:,1),nonLympCentroids(:,2));
            % Save the masks
            imwrite(lymp_mask,[mask_lymp_path patch_name]);
            imwrite(nonlymp_mask,[mask_nonlymp_path patch_name]);
        end
    end
    %disp('calculating the feature')
    try
        % Get the basics from the image and masks
        [lympCentroids,lympnucFeatures] = getNucLocalSimpleFeatures(patch_image,lymp_mask);
        [nonLympCentroids,~] = getNucLocalSimpleFeatures(patch_image,nonlymp_mask);
        % separate the nuclei area feature
        lympnucAreas=lympnucFeatures(:,1);
        % Get the dentil feature
        %[denFeat]=getDenTILFeatures(patch_image,lympCentroids,nonLympCentroids,lympnucAreas);
        coords={lympCentroids,nonLympCentroids,};
        [spaFeat_2gr_37,~]=getSpaTILFeatures_v2(coords);
        [spaFeat_2gr_40,~]=getSpaTILFeatures_v2(coords);
        [spaFeat_2gr_42,~]=getSpaTILFeatures_v2(coords);
        [spaFeat_2gr_43,~]=getSpaTILFeatures_v2(coords);
        % combine features
        spaFeat_2gr_alphas = [spaFeat_2gr_37 spaFeat_2gr_40 spaFeat_2gr_42 spaFeat_2gr_43];
        % Save the denTIL features
        save([spatil_feature_path patch_name(1:end-4) '.mat'],'spaFeat_2gr_alphas')
    catch
        %disp('Feature could not be extracted')
    end
    pb.print(x,length(img_dir))
end
fprintf('All done')
fprintf('\n')
