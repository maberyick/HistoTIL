clear all
close all
debug=0;

%%
% %personal computer paths
% % dbstop if error
originatrimfilepath='F:\Jon\Breast Cancer Data\S2_diced\S2_bulk\WideFOV\originals\';
% maskpath='F:\Jon\Breast Cancer Data\S2_diced\DL_nucleimasks\hhll_clean_oct14\';
% maskpath='F:\Jon\Breast Cancer Data\S2_diced\DL_nucleimasks\erosion_approach_masks\';
maskpath='F:\Jon\Breast Cancer Data\S2_diced\S2_bulk\WideFOV\nuclear_masks\';
epstromaskpath='F:\Jon\Breast Cancer Data\S2_diced\S2_bulk\WideFOV\epstroma_masks\';
% xylist='F:\Jon\Breast Cancer Data\S2_diced\data\8863-xy named';
exceldocpath='F:\Jon\Breast Cancer Data\Oncotype scores for S2 for r.xlsx';
load('F:\Jon\Jons Repository Code\Histology_models\Lymphocyte_model_german.mat');
% 
outputpath='F:\Jon\Breast Cancer Data\Feature output\WideFOV_July_2016_3.mat';
% 
boundary_list_not_mask=0;

if boundary_list_not_mask==1
    load('F:\Jon\Breast Cancer Data\Feature output\batch_bounds.mat');
end

pathsave='F:\Jon\Breast Cancer Data\Feature output\temp.mat';

%%
% HPC paths
% 
% originatrimfilepath='/mnt/pan/Data7/jrw178/S2_data/S2_big/originals/';
% % maskpath='F:\Jon\Breast Cancer Data\S2_diced\DL_nucleimasks\hhll_clean_oct14\';
% % maskpath='F:\Jon\Breast Cancer Data\S2_diced\DL_nucleimasks\erosion_approach_masks\';
% maskpath='/mnt/pan/Data7/jrw178/S2_data/S2_big/nuclear_mask/';
% epstromaskpath='/mnt/pan/Data7/jrw178/S2_data/S2_big/epstroma_masks/';
% % xylist='F:\Jon\Breast Cancer Data\S2_diced\data\8863-xy named';
% exceldocpath='/mnt/pan/Data7/jrw178/excel/Oncotype scores for S2 for r.xlsx';
% outputpath='/mnt/pan/Data7/jrw178/S2_data/S2_big/WideFOV_July_2016.mat';
% load('/mnt/Data7/jrw178/code/Lymphocyte_model_german.mat');
% pathsave='/mnt/pan/Data7/jrw178/S2_data/S2_big/temp.mat';

boundary_list_not_mask=0;

if boundary_list_not_mask==1
    load('F:\Jon\Breast Cancer Data\Feature output\batch_bounds.mat');
end

%% 
if exist(pathsave, 'file')
    load(pathsave);
end
if ~exist('trimfilepath', 'var')

originals = dir(originatrimfilepath);
if boundary_list_not_mask==0
    maskpaths = dir(maskpath);
end
% load(xylist)
[ndata, text, alldata] = xlsread(exceldocpath);
% [lhdata, ~, ~] = xlsread(llhhdocpath);
j=1;

%% This section gets a list of original files that have corresponding nuclei masks
% matlabpool open 8
folder1 = dir([originatrimfilepath,'*.png']);
name1 = {folder1.name};
for i=1:length(name1)
    filelist(i).paths=strcat(originatrimfilepath, name1{i});
    [~,name,~]=fileparts(filelist(i).paths);
    filelist(i).name=name;
end
[epstrolist]=Get_files_with_extension_func(epstromaskpath, '.png');
% for i=1:length(originals)
%     if ~isempty(strfind(originals(i).name, 'tif'));
%         filepath(j).original=strcat(originatrimfilepath, originals(i).name);
%         [~,filepath(j).originalname,~]=fileparts(filepath(j).original);
%         j=j+1;
%     end
% end
for i=1:length(filelist)
    if boundary_list_not_mask==0
        filelist(i).mask=strcat([maskpath, filelist(i).name, '_mask.png']);
    end
    filelist(i).esprob=[epstromaskpath, filelist(i).name, '.png'];
end
% for i=1:length(filelist)
%     for j=1:length(maskpaths)
%         if ~isempty(findstr(maskpaths(j).name, 'mask'));
%             maskname=strsplit(maskpaths(j).name, '_mask');
%             if ~isempty(strfind(maskpaths(j).name, filelist(i).name));
%                 filelist(i).mask=strcat(maskpath, maskpaths(j).name);
%                 filelist(i).prob=strcat(maskpath, maskname{1}, '_class.png');
%                 filelist(i).esprob=[epstromaskpath, maskname{1}, '_class.png'];
%             end
%         end
%     end
% 
% end

j=1; erchk=0;
for i=1:length(filelist)
    if ~isempty(filelist(i).mask);
        mimg=imread(filelist(i).mask);
        for k=1:length(ndata(:,1))
            if (max(max(mimg))>0 && ~isempty(findstr(num2str(ndata(k,1)), filelist(i).name)))
                if exist(filelist(i).paths, 'file')
                    trimfilepath(j).original=filelist(i).paths;
                else
                    disp(['Missing ', filelist(i).paths])
                    dbstop
                end
                trimfilepath(j).originalname=filelist(i).name;
                if boundary_list_not_mask==0
                    if exist(filelist(i).mask, 'file')
                        trimfilepath(j).mask=filelist(i).mask;
                    else
                        disp(['Missing ', filelist(i).mask])
                        dbstop
                    end
                end
                trimfilepath(j).odx=ndata(k,23);
%                 trimfilepath(j).prob=filelist(i).prob;
                if exist(filelist(i).esprob, 'file')
                    trimfilepath(j).esprob=filelist(i).esprob;
                else
                    disp(['Missing ', filelist(i).esprob])
                    dbstop
                end
                j=j+1;
            end
        end
    end
end

%% check to see which ones trimfilepath didn't get

k3=0;

for k1=1:length(filelist)
    m=0;
    for k2=1:length(trimfilepath)
        if strcmp(trimfilepath(k2).originalname, filelist(k1).name)
            m=1;
        end
    end
    if m==0
        k3=k3+1;
        missinglist(k3).name=filelist(k1).name;
    end
end
   save(pathsave, 'trimfilepath');     
end
% [namelist]=get_unique_cell_vector(epdatabase(2:end,end-1));

disp('paths set up');
%% This section gets the bounds and texture features
         
% bigbounds.r=[]; bigbounds.c=[]; bigbounds.centroid_r=[]; bigbounds.centroid_c=[];

% [~, description]=texturefeats(trimfilepath(1).original, trimfilepath(1).mask);
% texturesets=cell(length(trimfilepath), length(description));

%% This should check to see if output already exists, and attempt to pick up where the last run left off if it does

if (exist(outputpath, 'file') ~= 0)
    load(outputpath);
    [featureheight, featurewidth]=size(epdatabase);
    [namelist]=get_unique_cell_vector(epdatabase(2:end,end-1));
    filestart = [];
    
%% use for missing datapoints
if (exist ('epdatabase', 'var') ~= 0)
    for i=1:length(trimfilepath)
        if ~(ismember(trimfilepath(i).originalname, namelist))
            filestart=[filestart, i];
        end
    end
else
    filestart=2;
end 
    
%% use for parallel processing    
%     if (exist ('epdatabase', 'var') ~= 0)
%         [featureheight, featurewidth]=size(epdatabase);
%         filestart = [];
%         for i=1:featureheight-2
%             if isempty(epdatabase{i+1,1})
%                 filestart=[filestart,i+1];
%             end
%         end
%     else
%         filestart=2:featureheight;
%     end
%%
else %Run a test run to find out how many features are involved
%     lowdatabase=cell({});
    [img, epitheliummask, stromamask]=paths2esboundsnimg_2(trimfilepath, 1);  %CAUTION EXTREME FILTERING May 2016
%     [epitheliummask, stromamask, ep_lymphmask, st_lymphmask]=separate_lymphocyte_masks_2(img, epitheliummask, stromamask, model);
    [tempepfeats, tempstrofeats, labels]=epstroma_features_2(img, epitheliummask, stromamask, trimfilepath, 1);  %CAUTION: architecture feats turned off May 2016
%     [tempepfeats, tempstrofeats]=add_lymph_features(tempepfeats, tempstrofeats, ep_lymphmask, st_lymphmask, epitheliummask, stromamask);
     featureheight=length(trimfilepath)+1;
    featurewidth=size(tempepfeats,2);
    epdatabase=cell(featureheight, featurewidth);
    stdatabase=cell(featureheight, featurewidth);
    for k1=1:featurewidth
        epdatabase{1,k1}=tempepfeats{1,k1};
        stdatabase{1,k1}=tempepfeats{1,k1};
    end
    filestart=2:featureheight;
    save(outputpath, 'epdatabase', 'stdatabase');
end
disp('Labels have been initialized');

%% start the feature extraction

%% get parallel list set up
psz=4;
for k1=1:floor(length(filestart)/psz)
    parset(k1).list=filestart((k1-1)*psz+1:k1*psz);  
end
if exist('parset','var')
    parset(end+1).list=filestart(k1*psz+1:end);
else
    parset.list=[];
end

%%

if ~(isempty(filestart))
    disp(strcat('Resuming work at line _', num2str(filestart(1))));
   %Getting the i is tricky.  We're going to say that featureheight is the
   %number of samples +1 for the labels.
   %That means that i needs to be replaced with i-1 when referring to the
   %trimfilepath number.
bigtic=tic;
 for k1=1:length(parset)
       lowtic=tic;
       try
        parfor i=parset(k1).list % This function can be a parfor loop if you comment out the save function
        disp(['running_', trimfilepath(i-1).originalname]);
            [img, epitheliummask, stromamask]=paths2esboundsnimg_2(trimfilepath, i-1);
%             [epitheliummask, stromamask, ep_lymphmask, st_lymphmask]=separate_lymphocyte_masks_2(img, epitheliummask, stromamask, model);
            [epcellout, stcellout, labeltemp]=epstroma_features_2(img, epitheliummask, stromamask, trimfilepath, i-1);
%             [epcellout, stcellout]=add_lymph_features(epcellout, stcellout, ep_lymphmask, st_lymphmask, epitheliummask, stromamask);
            for k2=1:featurewidth
                epdatabase{i,k2}=epcellout{2,k2};
                stdatabase{i,k2}=stcellout{2,k2};
            end
            
        end
         disp([num2str(toc(bigtic)/60), ' min in.  ', num2str((length(parset)-k1)*psz), ' images left.']);
        minperpic=(toc(lowtic)/60)/(psz);
        disp([num2str(minperpic), ' min/pic']);
        minleft=(length(parset)-k1)*psz*minperpic;
        if minleft > 60
            disp([num2str(minleft/60), ' hours remaining']);
        else
            disp([num2str(minleft), ' min remaining']);
        end
        
        save(outputpath, 'epdatabase', 'stdatabase');
       catch
        for i=parset(k1).list % This function can be a parfor loop if you comment out the save function
        disp(['running_', trimfilepath(i-1).originalname]);
            [img, epitheliummask, stromamask]=paths2esboundsnimg_2(trimfilepath, i-1);
%             [epitheliummask, stromamask, ep_lymphmask, st_lymphmask]=separate_lymphocyte_masks_2(img, epitheliummask, stromamask, model);
            [epcellout, stcellout, labeltemp]=epstroma_features_2(img, epitheliummask, stromamask, trimfilepath, i-1);
%             [epcellout, stcellout]=add_lymph_features(epcellout, stcellout, ep_lymphmask, st_lymphmask, epitheliummask, stromamask);
            for k2=1:featurewidth
                epdatabase{i,k2}=epcellout{2,k2};
                stdatabase{i,k2}=stcellout{2,k2};
            end    
        end
         disp([num2str(toc(bigtic)/60), ' min in.  ', num2str((length(parset)-k1)*psz), ' images left.']);
        minperpic=(toc(lowtic)/60)/(psz);
        disp([num2str(minperpic), ' min/pic']);
        minleft=(length(parset)-k1)*psz*minperpic;
        if minleft > 60
            disp([num2str(minleft/60), ' hours remaining']);
        else
            disp([num2str(minleft), ' min remaining']);
        end
        save(outputpath, 'epdatabase', 'stdatabase');
       end
   end
end
disp('Feature extraction complete');



% labeltemp=labels;
% for i=1:length(labeltemp)
%     epdatabase{1,i}=labeltemp{i};
%     stdatabase{1,i}=labeltemp{i};
% end

% disp('Labels added to final databases');
%% This computes large scale features

% [featslargescale,descriptionlargescale] = georges_largescale_features(bigbounds);
% outparameters=featurecell2array(featslargescale,descriptionlargescale);
% 
% [~,temp]=size(outparameters);
% for i=1:length(texturefeats)
%     outparameters{1,temp+i}=description{1,1}{1,i};
%     outparameters{2,temp+i}=texturefeats(i);
% end

%% save the output
[h1,~]=size(epdatabase);
[h2,~]=size(stdatabase);


for k1=1:featurewidth
    for k2=1:h1
        totaldata{k2,k1}=epdatabase{k2,k1};
    end
    top=h1;
    for k2=2:h2
        totaldata{k2+top-1,k1}=stdatabase{k2,k1};
    end
end

save(outputpath, 'epdatabase', 'stdatabase', 'totaldata');
%[feature_concat]=compute_CoLlAGe2D(img(:,:,1), epitheliummask, 2);
% matlabpool close
disp('Done!');