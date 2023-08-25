% Gather the features and ave for each cohort
%% add path of the dependencies
% Windows
addpath(genpath('C:\Users\cbarr23\Documents\HistoTIL'))
% Linux
%addpath(genpath('/home/maberyick/CCIPD_Research/Github/HistoTIL'))
%% paths and names
cohort_name = 'ccf';
% Windows
% E:\sclc\data\uh
% D:\Projects\HistoTIL\histotil_features\upenn
folder_matpatches = ['D:\Projects\HistoTIL\histotil_features\' cohort_name '\dataset_output\'];
save_path_full = ['D:\Projects\HistoTIL\histotil_features\' cohort_name '\histotil_features\'];
% Linux
%folder_matpatches = ['/media/maberyick/Elements/CCIPD_Projects/PhenoTIL_IO/' cohort_name '/histotil/histotil_features/dataset_output/'];
%save_path_full = ['/media/maberyick/Elements/CCIPD_Projects/PhenoTIL_IO/' cohort_name '/histotil/histotil_features/'];
mkdir(save_path_full)
folderList = dir(folder_matpatches);
folderNames = {folderList([folderList.isdir]).name};
folderNames = folderNames(~ismember(folderNames ,{'.','..'}));
folderNames_count = folderNames;
nuc_all_cohort_feat = [];
nuc_nontil_cohort_feat = [];
cot_all_cohort_feat = [];
cot_nontil_cohort_feat = [];
den_cohort_feat = [];
spa_cohort_feat = [];
spd_cohort_feat = [];
per_cohort_feat = [];

% Patch based
nuc_all_cohort_feat_patch = [];
nuc_nontil_cohort_feat_patch = [];
cot_all_cohort_feat_patch = [];
cot_nontil_cohort_feat_patch = [];
den_cohort_feat_patch = [];
spa_cohort_feat_patch = [];
spd_cohort_feat_patch = [];
per_cohort_feat_patch = [];
folderName_clean_patch = {};
counter_patch = 1;
%
folderName_clean = {};
counter = 1;

% loop through the values
textprogressbar('calculating feature stats: ');

%% Gather the local cell features
for mm=1:length(folderNames_count)
    perct = 100*mm/length(folderNames_count);
    textprogressbar(perct);
    nuc_all_feat_cohort = [];
    nuc_nontil_feat_cohort = [];
    contx_all_feat_cohort = [];
    contx_nontil_feat_cohort = [];
    folderName = folderNames_count{mm};
    path_comp = [folder_matpatches folderName filesep 'TIL_features' filesep];
    sub_folderList = dir([path_comp '*.mat']);
    if isempty(sub_folderList)
        folderNames(mm) = [];
        continue
    end
    for nn=1:length(sub_folderList)
        histotil_tmp = load([path_comp sub_folderList(nn).name]);
%         %% gather Nuclei features
        [nuc_all_feat_set, nuc_all_feat_name] = get_feature_set( ...
            histotil_tmp.nucFeatures, ...
            histotil_tmp.nucFeatures_epi, ...
            histotil_tmp.nucFeatures_stro, ...
            histotil_tmp.nucFeatures_bund, ...
            'histotil','nuclei',{'tiss','epi','stro','bund'});
        nuc_all_feat_cohort = [nuc_all_feat_cohort; nuc_all_feat_set];
        %% gather Nuclei features - non-TILs
        try
            histotil_tmp_nuc_feat_epi_nontil = histotil_tmp.nucFeatures_epi(~histotil_tmp.isLymphocyte_epi,:);
            if isempty(histotil_tmp_nuc_feat_epi_nontil)
                histotil_tmp_nuc_feat_epi_nontil = zeros(5,100);
            end
        catch
            histotil_tmp_nuc_feat_epi_nontil = zeros(5,100);
        end
        try
            histotil_tmp_nuc_feat_stro_nontil = histotil_tmp.nucFeatures_stro(~histotil_tmp.isLymphocyte_stro,:);
            if isempty(histotil_tmp_nuc_feat_stro_nontil)
                histotil_tmp_nuc_feat_stro_nontil = zeros(5,100);
            end
        catch
            histotil_tmp_nuc_feat_stro_nontil = zeros(5,100);
        end
        try
            histotil_tmp_nuc_feat_bund_nontil = histotil_tmp.nucFeatures_bund(~histotil_tmp.isLymphocyte_bund,:);
            if isempty(histotil_tmp_nuc_feat_bund_nontil)
                histotil_tmp_nuc_feat_bund_nontil = zeros(5,100);
            end
        catch
            histotil_tmp_nuc_feat_bund_nontil = zeros(5,100);
        end
        [nuc_nontil_feat_set, nuc_nontil_feat_name] = get_feature_set( ...
            histotil_tmp.nucFeatures(~histotil_tmp.isLymphocyte,:), ...
            histotil_tmp_nuc_feat_epi_nontil, ...
            histotil_tmp_nuc_feat_stro_nontil, ...
            histotil_tmp_nuc_feat_bund_nontil, ...
            'histotil','nuclei',{'tiss','epi','stro','bund'});
        nuc_nontil_feat_cohort = [nuc_nontil_feat_cohort; nuc_nontil_feat_set];
        %% gather contextual features
        [contx_all_feat_set, contx_all_feat_name] = get_feature_set( ...
            histotil_tmp.ctxFeat, ...
            histotil_tmp.ctxFeat_epi, ...
            histotil_tmp.ctxFeat_stro, ...
            histotil_tmp.ctxFeat_bund, ...
            'histotil','context',{'tiss','epi','stro','bund'});
        contx_all_feat_cohort = [contx_all_feat_cohort; contx_all_feat_set];
        % gather contextual features - non-TILs
        try
            histotil_tmp_contx_feat_epi_nontil = histotil_tmp.ctxFeat_epi(~histotil_tmp.isLymphocyte_epi,:);
            if isempty(histotil_tmp_contx_feat_epi_nontil)
                histotil_tmp_contx_feat_epi_nontil = zeros(5,87);
            end
        catch
            histotil_tmp_contx_feat_epi_nontil = zeros(5,87);
        end
        try
            histotil_tmp_contx_feat_stro_nontil = histotil_tmp.ctxFeat_stro(~histotil_tmp.isLymphocyte_stro,:);
            if isempty(histotil_tmp_contx_feat_stro_nontil)
                histotil_tmp_contx_feat_stro_nontil = zeros(5,87);
            end
        catch
            histotil_tmp_contx_feat_stro_nontil = zeros(5,87);
        end
        try
            histotil_tmp_contx_feat_bund_nontil = histotil_tmp.ctxFeat_bund(~histotil_tmp.isLymphocyte_bund,:);
            if isempty(histotil_tmp_contx_feat_bund_nontil)
                histotil_tmp_contx_feat_bund_nontil = zeros(5,87);
            end
        catch
            histotil_tmp_contx_feat_bund_nontil = zeros(5,87);
        end
        [contx_nontil_feat_set, contx_nontil_feat_name] = get_feature_set( ...
            histotil_tmp.ctxFeat(~histotil_tmp.isLymphocyte,:), ...
            histotil_tmp_contx_feat_epi_nontil, ...
            histotil_tmp_contx_feat_stro_nontil, ...
            histotil_tmp_contx_feat_bund_nontil, ...
            'histotil','context',{'tiss','epi','stro','bund'});
        contx_nontil_feat_cohort = [contx_nontil_feat_cohort; contx_nontil_feat_set];
    end
    %% gather again the local cell features to a higher level
    nuc_all_feat_cohort_total = [];
    nuc_nontil_feat_cohort_total = [];
    contx_all_feat_cohort_total = [];
    contx_nontil_feat_cohort_total = [];
    [siz_a_nuc,siz_b_nuc] = size(nuc_all_feat_cohort);
    [siz_a_nuc_nontil,siz_b_nuc_nontil] = size(nuc_nontil_feat_cohort);
    [siz_a_ctx,siz_b_ctx] = size(contx_nontil_feat_cohort);
    [siz_a_ctx_nontil,siz_b_ctx_nontil] = size(contx_nontil_feat_cohort);
    stats_size = 18;
    % gather patch based
    nuc_all_cohort_feat_patch = [nuc_all_cohort_feat_patch; nuc_all_feat_cohort];
    nuc_nontil_cohort_feat_patch = [nuc_nontil_cohort_feat_patch; nuc_nontil_feat_cohort];
    cot_all_cohort_feat_patch = [cot_all_cohort_feat_patch; contx_all_feat_cohort];
    cot_nontil_cohort_feat_patch = [cot_nontil_cohort_feat_patch; contx_nontil_feat_cohort];
    % ---------------------------------------------------- %
    if siz_a_nuc == 1
        for p=1:stats_size
            nuc_all_feat_cohort_total =[nuc_all_feat_cohort_total nuc_all_feat_cohort];
        end
    % if more than one feature point, calculate stats
    else
        nuc_all_feat_cohort_total = stats_calc_feature(nuc_all_feat_cohort);
    end
    % ---------------------------------------------------- %
    if siz_a_nuc_nontil == 1
        for p=1:stats_size
            nuc_nontil_feat_cohort_total =[nuc_nontil_feat_cohort_total nuc_nontil_feat_cohort];
        end
    % if more than one feature point, calculate stats
    else
        nuc_nontil_feat_cohort_total = stats_calc_feature(nuc_nontil_feat_cohort);
    end
    % ---------------------------------------------------- %
    if siz_a_ctx == 1
        for p=1:stats_size
            contx_all_feat_cohort_total =[contx_all_feat_cohort_total contx_all_feat_cohort];
            contx_nontil_feat_cohort_total =[contx_nontil_feat_cohort_total contx_nontil_feat_cohort];
        end
        % if more than one feature point, calculate stats
    else
        contx_all_feat_cohort_total = stats_calc_feature(contx_all_feat_cohort);
        contx_nontil_feat_cohort_total = stats_calc_feature(contx_nontil_feat_cohort);
    end
    % ---------------------------------------------------- %
    if siz_a_ctx_nontil == 1
        for p=1:stats_size
            contx_nontil_feat_cohort_total =[contx_nontil_feat_cohort_total contx_nontil_feat_cohort];
        end
        % if more than one feature point, calculate stats
    else
        contx_nontil_feat_cohort_total = stats_calc_feature(contx_nontil_feat_cohort);
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
    %
    spat2_feat_cohort_tiss = [];
    spat2_feat_cohort_epi = [];
    spat2_feat_cohort_stro = [];
    spat2_feat_cohort_bund = [];
    %
    perin_feat_cohort_tiss = [];
    perin_feat_cohort_epi = [];
    perin_feat_cohort_stro = [];
    perin_feat_cohort_bund = [];
    %
    for nn=1:length(sub_folderList)
        histotil_tmp = load([path_comp sub_folderList(nn).name]);
        %% gather density features
        if length(histotil_tmp.denFeat) == 19
            dens_feat_cohort_tiss = [dens_feat_cohort_tiss; histotil_tmp.denFeat];
        else
            dens_feat_cohort_tiss = [dens_feat_cohort_tiss; zeros(1,19)];
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
        %% gather spatial features v2
        if length(histotil_tmp.spaFeat_v2) == 1400
            spat2_feat_cohort_tiss = [spat2_feat_cohort_tiss; histotil_tmp.spaFeat_v2];
        else
            spat2_feat_cohort_tiss = [spat2_feat_cohort_tiss; zeros(1,1400)];
        end
        %
        if length(histotil_tmp.spaFeat_v2_epi) == 1400
            spat2_feat_cohort_epi = [spat2_feat_cohort_epi; histotil_tmp.spaFeat_v2_epi];
        else
            spat2_feat_cohort_epi = [spat2_feat_cohort_epi; zeros(1,1400)];
        end
        %
        if length(histotil_tmp.spaFeat_v2_stro) == 1400
            spat2_feat_cohort_stro = [spat2_feat_cohort_stro; histotil_tmp.spaFeat_v2_stro];
        else
            spat2_feat_cohort_stro = [spat2_feat_cohort_stro; zeros(1,1400)];
        end
        %
        if length(histotil_tmp.spaFeat_v2_bund) == 1400
            spat2_feat_cohort_bund = [spat2_feat_cohort_bund; histotil_tmp.spaFeat_v2_bund];
        else
            spat2_feat_cohort_bund = [spat2_feat_cohort_bund; zeros(1,1400)];
        end
        %
        %% gather peri-nuclear features
        if length(histotil_tmp.spaFeat_v2) == 1400
            perin_feat_cohort_tiss = [perin_feat_cohort_tiss; histotil_tmp.perinuclear];
        else
            perin_feat_cohort_tiss = [perin_feat_cohort_tiss; zeros(1,1400)];
        end
        %
        if length(histotil_tmp.spaFeat_v2_epi) == 1400
            perin_feat_cohort_epi = [perin_feat_cohort_epi; histotil_tmp.perinuclear_epi];
        else
            perin_feat_cohort_epi = [perin_feat_cohort_epi; zeros(1,1400)];
        end
        %
        if length(histotil_tmp.spaFeat_v2_stro) == 1400
            perin_feat_cohort_stro = [perin_feat_cohort_stro; histotil_tmp.perinuclear_stro];
        else
            perin_feat_cohort_stro = [perin_feat_cohort_stro; zeros(1,1400)];
        end
        %
        if length(histotil_tmp.spaFeat_v2_bund) == 1400
            perin_feat_cohort_bund = [perin_feat_cohort_bund; histotil_tmp.perinuclear_bund];
        else
            perin_feat_cohort_bund = [perin_feat_cohort_bund; zeros(1,1400)];
        end
        %
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
        'histotil','dentil',{'tiss','epi','stro','bund'});
    %% gather Spatial features
    [spat_feat_set_all, spat_feat_name] = get_feature_set( ...
        spat_feat_cohort_tiss, ...
        spat_feat_cohort_epi, ...
        spat_feat_cohort_stro, ...
        spat_feat_cohort_bund, ...
        'histotil','spatil',{'tiss','epi','stro','bund'});
    %% gather Spatial v2 features
    [spat2_feat_set_all, spat2_feat_name] = get_feature_set( ...
        spat2_feat_cohort_tiss, ...
        spat2_feat_cohort_epi, ...
        spat2_feat_cohort_stro, ...
        spat2_feat_cohort_bund, ...
        'histotil','spatil2',{'tiss','epi','stro','bund'});
    %% gather peri-nuclear features
    [peri_feat_set_all, peri_feat_name] = get_feature_set( ...
        perin_feat_cohort_tiss, ...
        perin_feat_cohort_epi, ...
        perin_feat_cohort_stro, ...
        perin_feat_cohort_bund, ...
        'histotil','perinuc',{'tiss','epi','stro','bund'});
    % Gather density and spatial Patch based
    den_cohort_feat_patch = [den_cohort_feat_patch; [dens_feat_cohort_tiss dens_feat_cohort_epi dens_feat_cohort_stro dens_feat_cohort_bund]];
    spa_cohort_feat_patch = [spa_cohort_feat_patch; [spat_feat_cohort_tiss spat_feat_cohort_epi spat_feat_cohort_stro spat_feat_cohort_bund]];
    spd_cohort_feat_patch = [spd_cohort_feat_patch; [spat2_feat_cohort_tiss spat2_feat_cohort_epi spat2_feat_cohort_stro spat2_feat_cohort_bund]];
    per_cohort_feat_patch = [per_cohort_feat_patch; [perin_feat_cohort_tiss perin_feat_cohort_epi perin_feat_cohort_stro perin_feat_cohort_bund]];
    % Gather the features for the cohort
    nuc_all_cohort_feat = [nuc_all_cohort_feat; nuc_all_feat_cohort_total];
    nuc_nontil_cohort_feat = [nuc_nontil_cohort_feat; nuc_nontil_feat_cohort_total];
    cot_all_cohort_feat = [cot_all_cohort_feat; contx_all_feat_cohort_total];
    cot_nontil_cohort_feat = [cot_nontil_cohort_feat; contx_nontil_feat_cohort_total];
    den_cohort_feat = [den_cohort_feat; dens_feat_cohort_all];
    spa_cohort_feat = [spa_cohort_feat; spat_feat_set_all];
    spd_cohort_feat = [spd_cohort_feat; spat2_feat_set_all];
    per_cohort_feat = [per_cohort_feat; peri_feat_set_all];
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
len_feat_name = length(nuc_nontil_feat_name);
total_feats_size = stats_size*len_feat_name;
varNames_local = strings(1,total_feats_size);
counter = 0;
for m=1:stats_size
    for k=1:length(nuc_nontil_feat_name)
        counter=counter+1;
        varNames_local{counter} = strcat(nuc_nontil_feat_name{k},'_var_',num2str(k),'_',varTypes{m});
    end
end

% contextual
len_feat_name = length(contx_nontil_feat_name);
total_feats_size = stats_size*len_feat_name;
varNames_contx = strings(1,total_feats_size);
counter = 0;
for m=1:stats_size
    for k=1:length(contx_nontil_feat_name)
        counter=counter+1;
        varNames_contx{counter} = strcat(contx_nontil_feat_name{k},'_var_',num2str(k),'_',varTypes{m});
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
len_feat_name = length(contx_all_feat_name);
varNames_contx_patch = strings(1,len_feat_name);
counter = 0;
for k=1:length(contx_all_feat_name)
    counter=counter+1;
    varNames_contx_patch{counter} = strcat('contx_feat_',num2str(k));
end
fid = fopen([save_path_full cohort_name '_patch_contextual_features_names.txt'], 'wt');
fprintf(fid, '%s\n', contx_all_feat_name);
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
len_feat_name = length(nuc_all_feat_name);
varNames_local_patch = strings(1,len_feat_name);
counter = 0;
for k=1:length(nuc_all_feat_name)
    counter=counter+1;
    varNames_local_patch{counter} = strcat('nucl_feat_',num2str(k));
end
fid = fopen([save_path_full cohort_name '_patch_nuclei_features_names.txt'], 'wt');
fprintf(fid, '%s\n', nuc_all_feat_name);
fclose(fid);
% density
varNames_density_inter = strings(1,length(dens_feat_name));
counter = 0;
for k=1:length(dens_feat_name)
    counter=counter+1;
    varNames_density_inter{counter} = strcat('density_feat_',num2str(k));
end
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
% spatial v2
varNames_spatial2_inter = strings(1,length(spat2_feat_name));
counter = 0;
for k=1:length(spat2_feat_name)
    counter=counter+1;
    varNames_spatial2_inter{counter} = strcat('spatial2_feat_',num2str(k));
end
%
fid = fopen([save_path_full cohort_name '_spatial2_features_names.txt'], 'wt');
fprintf(fid, '%s\n', spat2_feat_name);
fclose(fid);
% spatial v2 Patch
varNames_spatial2_inter_patch = strings(1,1400*4);
counter = 0;
for k=1:1400*4
    counter=counter+1;
    varNames_spatial2_inter_patch{counter} = strcat('spatial2_feat_',num2str(k));
end

% peri-nuclear
varNames_perinuc_inter = strings(1,length(peri_feat_name));
counter = 0;
for k=1:length(peri_feat_name)
    counter=counter+1;
    varNames_perinuc_inter{counter} = strcat('perinuc_feat_',num2str(k));
end
%
fid = fopen([save_path_full cohort_name '_perinuc_features_names.txt'], 'wt');
fprintf(fid, '%s\n', peri_feat_name);
fclose(fid);
% peri-nuclear Patch
varNames_perinuc_inter_patch = strings(1,37*4);
counter = 0;
for k=1:37*4
    counter=counter+1;
    varNames_perinuc_inter_patch{counter} = strcat('perinuc_feat_',num2str(k));
end
% Save the feature table
file_ID = folderNames';
% contextual
contextual_features = array2table(cot_all_cohort_feat,"VariableNames",varNames_contx_inter);
contextual_features = addvars(contextual_features,file_ID,'Before',1);
contextual_features.file_ID = string(contextual_features.file_ID);
% Save the tables
writetable(contextual_features,[save_path_full cohort_name '_contextual_features.csv']);

% contextual - nonTIL
contextual_features_nontil = array2table(cot_nontil_cohort_feat,"VariableNames",varNames_contx_inter);
contextual_features_nontil = addvars(contextual_features_nontil,file_ID,'Before',1);
contextual_features_nontil.file_ID = string(contextual_features_nontil.file_ID);
% Save the tables
writetable(contextual_features_nontil,[save_path_full cohort_name '_contextual_nontil_features.csv']);

% nuclei
nuclei_features = array2table(nuc_all_cohort_feat,"VariableNames",varNames_nucl_inter);
nuclei_features = addvars(nuclei_features,file_ID,'Before',1);
nuclei_features.file_ID = string(nuclei_features.file_ID);
% Save the tables
writetable(nuclei_features,[save_path_full cohort_name '_nuclei_features.csv']);

% nuclei - non-TIL
nuclei_features_nontil = array2table(nuc_nontil_cohort_feat,"VariableNames",varNames_nucl_inter);
nuclei_features_nontil = addvars(nuclei_features_nontil,file_ID,'Before',1);
nuclei_features_nontil.file_ID = string(nuclei_features_nontil.file_ID);
% Save the tables
writetable(nuclei_features_nontil,[save_path_full cohort_name '_nuclei_nontil_features.csv']);

% density
density_features = array2table(den_cohort_feat,"VariableNames",varNames_density_inter);
density_features = addvars(density_features,file_ID,'Before',1);
density_features.file_ID = string(density_features.file_ID);
% Save the tables
writetable(density_features,[save_path_full cohort_name '_density_features.csv']);

% spatial
spatial_features = array2table(spa_cohort_feat,"VariableNames",varNames_spatial_inter);
spatial_features = addvars(spatial_features,file_ID,'Before',1);
spatial_features.file_ID = string(spatial_features.file_ID);
% Save the tables
writetable(spatial_features,[save_path_full cohort_name '_spatial_features.csv']);

% spatial v2
spatial2_features = array2table(spd_cohort_feat,"VariableNames",varNames_spatial2_inter);
spatial2_features = addvars(spatial2_features,file_ID,'Before',1);
spatial2_features.file_ID = string(spatial2_features.file_ID);
% Save the tables
writetable(spatial2_features,[save_path_full cohort_name '_spatial2_features.csv']);

% peri-nuclear
pernuc_features = array2table(per_cohort_feat,"VariableNames",varNames_perinuc_inter);
pernuc_features = addvars(pernuc_features,file_ID,'Before',1);
pernuc_features.file_ID = string(pernuc_features.file_ID);
% Save the tables
writetable(pernuc_features,[save_path_full cohort_name '_perinuclear_features.csv']);