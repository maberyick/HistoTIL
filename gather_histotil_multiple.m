% Gather the features and ave for each cohort
%% add path of the dependencies
addpath(genpath('/home/maberyick/CCIPD_Research/Github/HistoTIL'))
%% paths and names
cohort_name = 'AZ';
folder_matpatches = ['/media/maberyick/Elements/CCIPD_Projects/PhenoTIL_IO/' cohort_name '/histotil/histotil_features/dataset_output/'];
save_path_full = ['/media/maberyick/Elements/CCIPD_Projects/PhenoTIL_IO/' cohort_name '/histotil/histotil_features/'];
folderList = dir(folder_matpatches);
folderNames = {folderList([folderList.isdir]).name};
folderNames = folderNames(~ismember(folderNames ,{'.','..'}));

nuc_cohort_feat = [];
cot_cohort_feat = [];
den_cohort_feat = [];
spa_cohort_feat = [];

% Patch based
nuc_cohort_feat_patch = [];
cot_cohort_feat_patch = [];
den_cohort_feat_patch = [];
spa_cohort_feat_patch = [];
folderName_clean_patch = {};
counter_patch = 1;
%
folderName_clean = {};
counter = 1;

% loop through the values
textprogressbar('calculating feature stats: ');

%% Gather the local cell features
for mm=1:length(folderNames)
    perct = 100*mm/length(folderNames);
    textprogressbar(perct);
    nuc_feat_cohort = [];
    contx_feat_cohort = [];
    folderName = folderNames{mm};
    path_comp = [folder_matpatches folderName '/TIL_features/'];
    sub_folderList = dir([path_comp '*.mat']);
    if isempty(sub_folderList)
        continue
    end
    for nn=1:length(sub_folderList)
        histotil_tmp = load([path_comp sub_folderList(nn).name]);
        %% gather Nuclei features
        [nuc_feat_set, nuc_feat_name] = get_feature_set( ...
            histotil_tmp.nucFeatures, ...
            histotil_tmp.nucFeatures_epi, ...
            histotil_tmp.nucFeatures_stro, ...
            histotil_tmp.nucFeatures_bund, ...
            'histotil','nuclei',{'tiss','epi','stro','bund'});
        nuc_feat_cohort = [nuc_feat_cohort; nuc_feat_set];
        %% gather contextual features
        [contx_feat_set, contx_feat_name] = get_feature_set( ...
            histotil_tmp.ctxFeat, ...
            histotil_tmp.ctxFeat_epi, ...
            histotil_tmp.ctxFeat_stro, ...
            histotil_tmp.ctxFeat_bund, ...
            'histotil','context',{'tiss','epi','stro','bund'});
        contx_feat_cohort = [contx_feat_cohort; contx_feat_set];    
    end
    %% gather again the local cell features to a higher level
    nuc_feat_cohort_total = [];
    contx_feat_cohort_total = [];
    [siz_a_nuc,siz_b_nuc] = size(nuc_feat_cohort);
    [siz_a_ctx,siz_b_ctx] = size(contx_feat_cohort);
    stats_size = 18;
    % gather patch based
    nuc_cohort_feat_patch = [nuc_cohort_feat_patch; nuc_feat_cohort];
    cot_cohort_feat_patch = [cot_cohort_feat_patch; contx_feat_cohort];
    % ---------------------------------------------------- %
    if siz_a_nuc == 1
        for p=1:stats_size
            nuc_feat_cohort_total =[nuc_feat_cohort_total nuc_feat_cohort];
        end
    % if more than one feature point, calculate stats
    else
        nuc_feat_cohort_total = stats_calc_feature(nuc_feat_cohort);
    end
    % ---------------------------------------------------- %
    if siz_a_ctx == 1
        for p=1:stats_size
            contx_feat_cohort_total =[contx_feat_cohort_total contx_feat_cohort];
        end
        % if more than one feature point, calculate stats
    else
        contx_feat_cohort_total = stats_calc_feature(contx_feat_cohort);
    end
    % ---------------------------------------------------- %
    %% gather the density and spatial features
    dens_feat_cohort_tiss = [];
    dens_feat_cohort_epi = [];
    dens_feat_cohort_stro = [];
    dens_feat_cohort_bund = [];
    %cot_cohort_feat_patch
    spat_feat_cohort_tiss = [];
    spat_feat_cohort_epi = [];
    spat_feat_cohort_stro = [];
    spat_feat_cohort_bund = [];
    for nn=1:length(sub_folderList)
        histotil_tmp = load([path_comp sub_folderList(nn).name]);
        %% gather density features
        if length(histotil_tmp.denFeat) == 19
            dens_feat_cohort_tiss = [dens_feat_cohort_tiss; histotil_tmp.denFeat];
        else
            dens_feat_cohort_tiss = [dens_feat_cohort_tiss; histotil_tmpzeros(1,19)];
        end
        %
        if length(histotil_tmp.denFeat_epi) == 19
            dens_feat_cohort_epi = [dens_feat_cohort_epi; histotil_tmp.denFeat_epi];
        else
            dens_feat_cohort_epi = [dens_feat_cohort_epi; zeros(1,19)];
        end
        %
        if length(histotil_tmp.denFeat_stro) == 19
            dens_feat_cohort_stro = [dens_feat_cohort_stro; histotil_tmp.denFeat_stro];
        else
            dens_feat_cohort_stro = [dens_feat_cohort_stro; zeros(1,19)];
        end
        %
        if length(histotil_tmp.denFeat_bund) == 19
            dens_feat_cohort_bund = [dens_feat_cohort_bund; histotil_tmp.denFeat_bund];
        else
            dens_feat_cohort_bund = [dens_feat_cohort_bund; zeros(1,19)];
        end
        %% gather spatial features
        if length(histotil_tmp.spaFeat) == 85
            spat_feat_cohort_tiss = [spat_feat_cohort_tiss; histotil_tmp.spaFeat];
        else
            spat_feat_cohort_tiss = [spat_feat_cohort_tiss; zeros(1,85)];
        end
        %
        if length(histotil_tmp.spaFeat_epi) == 85
            spat_feat_cohort_epi = [spat_feat_cohort_epi; histotil_tmp.spaFeat_epi];
        else
            spat_feat_cohort_epi = [spat_feat_cohort_epi; zeros(1,85)];
        end
        %
        if length(histotil_tmp.spaFeat_stro) == 85
            spat_feat_cohort_stro = [spat_feat_cohort_stro; histotil_tmp.spaFeat_stro];
        else
            spat_feat_cohort_stro = [spat_feat_cohort_stro; zeros(1,85)];
        end
        %
        if length(histotil_tmp.spaFeat_bund) == 85
            spat_feat_cohort_bund = [spat_feat_cohort_bund; histotil_tmp.spaFeat_bund];
        else
            spat_feat_cohort_bund = [spat_feat_cohort_bund; zeros(1,85)];
        end
        folderName_clean_patch{counter_patch} = folderName;
        counter_patch = counter_patch+1;
    end
    histotil_tmp = load([path_comp sub_folderList(nn).name]);
    %% gather Density features
    [dens_feat_cohort_all, dens_feat_name] = get_feature_set( ...
        dens_feat_cohort_tiss, ...
        dens_feat_cohort_epi, ...
        dens_feat_cohort_stro, ...
        dens_feat_cohort_bund, ...
        'histotil','nuclei',{'tiss','epi','stro','bund'});
    %% gather Spatial features
    [spat_feat_set_all, spat_feat_name] = get_feature_set( ...
        spat_feat_cohort_tiss, ...
        spat_feat_cohort_epi, ...
        spat_feat_cohort_stro, ...
        spat_feat_cohort_bund, ...
        'histotil','context',{'tiss','epi','stro','bund'});
    %% Gather density and spatial Patch based
    den_cohort_feat_patch = [den_cohort_feat_patch; [dens_feat_cohort_tiss dens_feat_cohort_epi dens_feat_cohort_stro dens_feat_cohort_bund]];
    spa_cohort_feat_patch = [spa_cohort_feat_patch; [spat_feat_cohort_tiss spat_feat_cohort_epi spat_feat_cohort_stro spat_feat_cohort_bund]];
    %% Gather the features for the cohort
    nuc_cohort_feat = [nuc_cohort_feat; nuc_feat_cohort_total];
    cot_cohort_feat = [cot_cohort_feat; contx_feat_cohort_total];
    den_cohort_feat = [den_cohort_feat; dens_feat_cohort_all];
    spa_cohort_feat = [spa_cohort_feat; spat_feat_set_all];
    folderName_clean{counter} = folderName;
    counter = counter+1;
end
fprintf('\n')
textprogressbar('done');
%% get the names for the local cell features
% nuclei
varTypes = {'stat_bds_1','stat_bds_2','stat_bds_3',...
    'stat_bds_4','stat_bds_5','stat_bds_6','stat_bds_7',...
    'stat_bds_8','stat_bds_9','stat_ctd_1','stat_ctd_2',...
    'stat_ctd_3','stat_ctd_4','stat_ctd_5','stat_ctd_6',...
    'stat_rdz_1','stat_rdz_2','stat_rdz_3'};
stats_size = length(varTypes);
len_feat_name = length(nuc_feat_name);
total_feats_size = stats_size*len_feat_name;
varNames_local = strings(1,total_feats_size);
counter = 0;
for m=1:stats_size
    for k=1:length(nuc_feat_name)
        counter=counter+1;
        varNames_local{counter} = strcat(nuc_feat_name{k},'_var_',num2str(k),'_',varTypes{m});
    end
end

% contextual
len_feat_name = length(contx_feat_name);
total_feats_size = stats_size*len_feat_name;
varNames_contx = strings(1,total_feats_size);
counter = 0;
for m=1:stats_size
    for k=1:length(contx_feat_name)
        counter=counter+1;
        varNames_contx{counter} = strcat(contx_feat_name{k},'_var_',num2str(k),'_',varTypes{m});
    end
end
%% Save the feature names and make them interpretable
% contextual
varNames_contx_inter = strings(1,length(varNames_contx));
counter = 0;
for k=1:length(varNames_contx)
    counter=counter+1;
    varNames_contx_inter{counter} = strcat('contx_feat_',num2str(k));
end
%
fid = fopen([save_path_full cohort_name '_contextual_features_names.txt'], 'wt');
fprintf(fid, '%s\n', varNames_contx);
fclose(fid);
% contextual Patch
len_feat_name = length(contx_feat_name);
varNames_contx_patch = strings(1,len_feat_name);
counter = 0;
for k=1:length(contx_feat_name)
    counter=counter+1;
    varNames_contx_patch{counter} = strcat('contx_feat_',num2str(k));
end
fid = fopen([save_path_full '/patch_based/' cohort_name '_patch_contextual_features_names.txt'], 'wt');
fprintf(fid, '%s\n', contx_feat_name);
fclose(fid);
% nuclei
varNames_nucl_inter = strings(1,length(varNames_local));
counter = 0;
for k=1:length(varNames_local)
    counter=counter+1;
    varNames_nucl_inter{counter} = strcat('nucl_feat_',num2str(k));
end
%
fid = fopen([save_path_full cohort_name '_nuclei_features_names.txt'], 'wt');
fprintf(fid, '%s\n', varNames_local);
fclose(fid);
% nuclei Patch
len_feat_name = length(nuc_feat_name);
varNames_local_patch = strings(1,len_feat_name);
counter = 0;
for k=1:length(nuc_feat_name)
    counter=counter+1;
    varNames_local_patch{counter} = strcat('nucl_feat_',num2str(k));
end
fid = fopen([save_path_full '/patch_based/' cohort_name '_patch_nuclei_features_names.txt'], 'wt');
fprintf(fid, '%s\n', nuc_feat_name);
fclose(fid);
% density
varNames_density_inter = strings(1,length(dens_feat_name));
counter = 0;
for k=1:length(dens_feat_name)
    counter=counter+1;
    varNames_density_inter{counter} = strcat('density_feat_',num2str(k));
end
fid = fopen([save_path_full cohort_name '_density_features_names.txt'], 'wt');
fprintf(fid, '%s\n', dens_feat_name);
fclose(fid);
% denstiy Patch
varNames_density_inter_patch = strings(1,19*4);
counter = 0;
for k=1:19*4
    counter=counter+1;
    varNames_density_inter_patch{counter} = strcat('density_feat_',num2str(k));
end
% spatial
varNames_spatial_inter = strings(1,length(spat_feat_name));
counter = 0;
for k=1:length(spat_feat_name)
    counter=counter+1;
    varNames_spatial_inter{counter} = strcat('spatial_feat_',num2str(k));
end
%
fid = fopen([save_path_full cohort_name '_spatial_features_names.txt'], 'wt');
fprintf(fid, '%s\n', spat_feat_name);
fclose(fid);
% spatial Patch
varNames_spatial_inter_patch = strings(1,85*4);
counter = 0;
for k=1:85*4
    counter=counter+1;
    varNames_spatial_inter_patch{counter} = strcat('spatial_feat_',num2str(k));
end
%% Save the feature table
% contextual
contextual_features = array2table(cot_cohort_feat,"VariableNames",varNames_contx_inter);
contextual_features_comp = addvars(contextual_features,folderName_clean','Before','contx_feat_1','NewVariableNames','case_ID');
% Save the tables
writetable(contextual_features_comp,[save_path_full cohort_name '_contextual_features.csv']);

% nuclei
nuclei_features = array2table(nuc_cohort_feat,"VariableNames",varNames_nucl_inter);
nuclei_features_comp = addvars(nuclei_features,folderName_clean','Before','nucl_feat_1','NewVariableNames','case_ID');
% Save the tables
writetable(nuclei_features_comp,[save_path_full cohort_name '_nuclei_features.csv']);

% density
density_features = array2table(den_cohort_feat,"VariableNames",varNames_density_inter);
density_features_comp = addvars(density_features,folderName_clean','Before','density_feat_1','NewVariableNames','case_ID');
% Save the tables
writetable(density_features_comp,[save_path_full cohort_name '_density_features.csv']);

% spatial
spatial_features = array2table(spa_cohort_feat,"VariableNames",varNames_spatial_inter);
spatial_features_comp = addvars(spatial_features,folderName_clean','Before','spatial_feat_1','NewVariableNames','case_ID');
% Save the tables
writetable(spatial_features_comp,[save_path_full cohort_name '_spatial_features.csv']);

%% Patch based
% contextual
contextual_features = array2table(cot_cohort_feat_patch,"VariableNames",varNames_contx_patch);
contextual_features_comp = addvars(contextual_features,folderName_clean_patch','Before','contx_feat_1','NewVariableNames','case_ID');
% Save the tables
writetable(contextual_features_comp,[save_path_full '/patch_based/' cohort_name '_contextual_features.csv']);

% nuclei
nuclei_features = array2table(nuc_cohort_feat_patch,"VariableNames",varNames_local_patch);
nuclei_features_comp = addvars(nuclei_features,folderName_clean_patch','Before','nucl_feat_1','NewVariableNames','case_ID');
% Save the tables
writetable(nuclei_features_comp,[save_path_full '/patch_based/' cohort_name '_nuclei_features.csv']);

% density
density_features = array2table(den_cohort_feat_patch,"VariableNames",varNames_density_inter_patch);
density_features_comp = addvars(density_features,folderName_clean_patch','Before','density_feat_1','NewVariableNames','case_ID');
% Save the tables
writetable(density_features_comp,[save_path_full '/patch_based/' cohort_name '_density_features.csv']);

% spatial
spatial_features = array2table(spa_cohort_feat_patch,"VariableNames",varNames_spatial_inter_patch);
spatial_features_comp = addvars(spatial_features,folderName_clean_patch','Before','spatial_feat_1','NewVariableNames','case_ID');
% Save the tables
writetable(spatial_features_comp,[save_path_full '/patch_based/' cohort_name '_spatial_features.csv']);