function [feat_cohort,varNames] = get_feature_set(feat_a,feat_b,feat_c,feat_d,feature_type,feature_name,feat_names)
%% get stats metrics from the feature and return a single row with all the features
%% Initiate the names of the stats and features
[feat_points_a,feat_initial_a] = size(feat_a);
[feat_points_b,~] = size(feat_b);
[feat_points_c,~] = size(feat_c);
[feat_points_d,~] = size(feat_d);
varTypes = {'stat_bds_1','stat_bds_2','stat_bds_3',...
            'stat_bds_4','stat_bds_5','stat_bds_6','stat_bds_7',...
            'stat_bds_8','stat_bds_9','stat_ctd_1','stat_ctd_2',...
            'stat_ctd_3','stat_ctd_4','stat_ctd_5','stat_ctd_6',...
            'stat_rdz_1','stat_rdz_2','stat_rdz_3'};
stats_size = length(varTypes);
len_feat_name = length(feat_names);
total_feats_size = feat_initial_a*stats_size*len_feat_name;
varNames = strings(1,total_feats_size);
counter = 0;
for n=1:len_feat_name
    for m=1:stats_size
        for k=1:feat_initial_a
            counter=counter+1;
            varNames{counter} = strcat(feature_type,'_',feature_name,'_',feat_names{n},'_var_',num2str(k),'_',varTypes{m});
        end
    end
end
%% Get the feature stats
% if single feature point, then reproduce across the feat columns
feat_file_group_a = [];
feat_file_group_b = [];
feat_file_group_c = [];
feat_file_group_d = [];
% ---------------------------------------------------- %
if feat_points_a == 1
    for p=1:stats_size
        feat_file_group_a =[feat_file_group_a feat_a];
    end
% if more than one feature point, calculate stats
else
    feat_file_group_a = stats_calc_feature(feat_a);
end
% ---------------------------------------------------- %
if feat_points_b == 1
    for p=1:stats_size
        feat_file_group_b =[feat_file_group_b feat_b];
    end
% if more than one feature point, calculate stats
else
    feat_file_group_b = stats_calc_feature(feat_b);
end
% ---------------------------------------------------- %
if feat_points_c == 1
    for p=1:stats_size
        feat_file_group_c =[feat_file_group_c feat_c];
    end
% if more than one feature point, calculate stats
else
    feat_file_group_c = stats_calc_feature(feat_c);
end
% ---------------------------------------------------- %
if feat_points_d == 1
    for p=1:stats_size
        feat_file_group_d =[feat_file_group_d feat_d];
    end
% if more than one feature point, calculate stats
else
    feat_file_group_d = stats_calc_feature(feat_d);
end
% ---------------------------------------------------- %

%% Return the variables
feat_cohort = [feat_file_group_a feat_file_group_b feat_file_group_c feat_file_group_d];
%%
end