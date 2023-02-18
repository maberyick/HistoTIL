% Gather the features and ave for each cohort
% add path of the dependencies
addpath(genpath('/home/maberyick/CCIPD_Research/Github/HistoTIL'))
cohort_name = 'AZ';
folder_matpatches = ['/media/maberyick/Elements/CCIPD_Projects/PhenoTIL_IO/' cohort_name '/histotil/histotil_features/dataset_output/'];
save_path_full = ['/media/maberyick/Elements/CCIPD_Projects/PhenoTIL_IO/' cohort_name '/histotil/histotil_features/'];
%save_path_full = '/home/maberyick/CCIPD_Research/Projects/pj_phenoTIL_Chemo_V1/phenoTIL_V1_rebuttal/dentil/';
%med_file_path = '/home/maberyick/CCIPD_Research/Projects/pj_phenoTIL_Chemo_V1/phenoTIL_V1_rebuttal/dentil/dentil_cohort_med_org.csv';
folderList = dir(folder_matpatches);
folderNames = {folderList([folderList.isdir]).name};
folderNames = folderNames(~ismember(folderNames ,{'.','..'}));


for mm=1:length(folderNames)
    folderName = folderNames{mm};
    path_comp = [folder_matpatches folderName '/TIL_features/'];
    sub_folderList = dir([path_comp '*.mat']);
    for nn=1:length(sub_folderList)
        histotil_tmp = load([path_comp sub_folderList(nn).name]);
        %% gather 
    end
end

%% Until here



feature_type = 'histotil';
cohort_dir = dir([cohort_path '*.mat']);
length_cohort = length(cohort_dir);
file_ID = repelem({'tmp'}, length_cohort)';
name_list_main_copy = repelem({'tmp'}, length_cohort)';
% Type of feature extracted from image (only for testing denTIL)
image_Type = "Single";
% if the features are already gathered per patient
feature_gather = "Yes";
% gather the original names of cases
if image_Type == "Multiple"
    for y=1:length(cohort_dir)
        file_name = cohort_dir(y).name;
        name_list_main_copy{y} = file_name;
        str_place = strfind(file_name,'_1.mat');
        file_main = file_name(1:str_place-1);
        file_ID{y} = file_main;
    end
elseif image_Type == "Single"
    for y=1:length(cohort_dir)
        file_name = cohort_dir(y).name;
        name_list_main_copy{y} = file_name;
        str_place = strfind(file_name,'.mat');
        file_main = file_name(1:str_place-1);
        file_ID{y} = file_main;
    end
end
file_ID = unique(file_ID);
% size of the initial feature set
% dentil
feat_initial = 19;
varTypes = {'stat_bds_1','stat_bds_2','stat_bds_3',...
            'stat_bds_4','stat_bds_5','stat_bds_6','stat_bds_7',...
            'stat_bds_8','stat_bds_9','stat_ctd_1','stat_ctd_2',...
            'stat_ctd_3','stat_ctd_4','stat_ctd_5','stat_ctd_6',...
            'stat_rdz_1','stat_rdz_2','stat_rdz_3'};
stats_size = length(varTypes);
total_feats_size = feat_initial*stats_size;
varNames = strings(1,total_feats_size);
counter = 0;
for m=1:stats_size
    for k=1:feat_initial
        counter=counter+1; 
        varNames{counter} = strcat( feature_type,'_var_',num2str(k),'_',varTypes{m} );
    end
end

length_cohort_main = length(file_ID);
feat_cohort = array2table([zeros(length_cohort_main,total_feats_size)],'VariableNames',varNames);
feat_cohort = addvars(feat_cohort,file_ID,'Before',varNames{1});

% loop through the values
textprogressbar('calculating feature stats: ');

for x=1:length_cohort_main
    perct = 100*x/length_cohort_main;
    textprogressbar(perct);
    file_name_main = file_ID{x};
    % find the files with the same name
    tmp_files_indx = startsWith(name_list_main_copy,file_name_main);
    tmp_files_pos = find(tmp_files_indx==1);
    file_comp = [];
    if feature_gather == "No"
        for y=1:length(tmp_files_pos)
            tmp_dentil_val = load([cohort_path name_list_main_copy{tmp_files_pos(y)}]);
            dentil_tmp = tmp_dentil_val.denFeat;
            file_comp = [file_comp; dentil_tmp];
        end
    elseif feature_gather == "Yes"
        for y=1:length(tmp_files_pos)
            tmp_dentil_val = load([cohort_path name_list_main_copy{tmp_files_pos(y)}]);
            % name given by German feature extraction (allDentilFeat)
            file_comp = tmp_dentil_val.allDentilFeat;
        end
    end
    % Get the mean value, as values are less than 25 and tend to work
    % better
    % if there is only one patch, no need to calculate stats.
    if length(tmp_files_pos) == 1
        if feature_gather == "No"
            feat_file_group = [];
            for p=1:stats_size
                feat_file_group =[feat_file_group file_comp];
            end
        elseif feature_gather == "Yes"
                feat_file_group = stats_calc_feature(file_comp);
        end
    else
        %dentil_file_mean = mean(dentil_file_comp);
        feat_file_group = stats_calc_feature(file_comp);
    end
    % save feature set to the table
    for k=1:total_feats_size
        eval(sprintf('feat_cohort.%s(%g) = feat_file_group(%g);',varNames{k},x,k));
    end
end
fprintf('\n')
textprogressbar('done');
% Save the tables
writetable(feat_cohort,[save_path cohort_name '.csv']);
%% Run this only once, if the file does not have the feature columns
% Read the med info excel file
%med_file = readtable(med_file_path);
% add empty columns of the variables
%emptyCol = cell(length(med_file.ID),length(varNames));
%emptyCol_name = array2table(emptyCol,'VariableNames',varNames);
%med_file_full = [med_file,emptyCol_name];
% Save the tables
%writetable(med_file_full,[save_path_full  'complete_list_v2.csv']);
%% save the features
med_file = readtable([save_path_full  'complete_list_v2.csv']);
feat_pos = find(strcmpi(med_file.Properties.VariableNames,varNames(1)));
for x=1:length_cohort_main
    try
        index = find(med_file.ID == string(feat_cohort.file_ID{x}));
        %med_file(index,:) = [med_file(index,feat_pos:end) feat_cohort(1,:)];
        med_file(index,feat_pos:end) = feat_cohort(x,2:end);
    catch
        disp('case file not found')
         string(feat_cohort.file_ID{x})
        continue
    end
end
writetable(med_file,[save_path_full  'complete_list_v2.csv']);
%or
% [~,index] = ismember(specific_value,TABLE{:,:})