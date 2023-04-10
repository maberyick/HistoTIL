%% paths and names
cohort_name = 'uh';
% Windows
%folder_matpatches = ['D:\Projects\PhenoTIL_IO\' cohort_name '\histotil\histotil_features\dataset_output\'];
%save_path_full = ['D:\Projects\PhenoTIL_IO\' cohort_name '\histotil\histotil_features\'];
% Linux
main_path = '/media/maberyick/TSPS2/sclc/data/';
folder_matpatches = [main_path cohort_name '/histotil_features/dataset_output/'];
folder_tilepatches = [main_path cohort_name '/tiles/'];
save_path_full = [main_path cohort_name '/histotil_features/'];
%
folderList = dir(folder_matpatches);
folderNames = {folderList([folderList.isdir]).name};
folderNames = folderNames(~ismember(folderNames ,{'.','..'}));
folderNames_count = folderNames;
%
z_val = 1.96;
p_val = 0.5;
e_val = 0.05;
b = (z_val^2)*(p_val*(1-p_val));
c = e_val^2;
a = b/c;
e = (z_val^2)*p_val*(1-p_val);
%
patches_cohort = [];
tiles_cohort = [];
tile_patch_per = [];
exclude_tiles = [];
tile_patch_per_2 = [];
neces_tiles = [];
folderName_clean_patch = {};
%% Gather the local cell features
for mm=1:length(folderNames_count)
    % count the mat files
    folderName = folderNames_count{mm};
    path_comp = [folder_matpatches folderName filesep 'TIL_features' filesep];
    sub_folderList = dir([path_comp '*.mat']);
    % count the tiles
    tiles_comp = [folder_tilepatches folderName filesep];
    sub_folderList_tile = dir([tiles_comp '*.png']);
    %
    folderName_clean_patch{mm} = folderName;
    patches_cohort = [patches_cohort; length(sub_folderList)];
    tiles_cohort = [tiles_cohort; length(sub_folderList_tile)];
    exclude_tiles = [exclude_tiles; length(sub_folderList_tile) - length(sub_folderList)];
    tile_patch_per = [tile_patch_per; 100*length(sub_folderList)/ length(sub_folderList_tile)];
    %
    % number of samples needed per image
    f = length(sub_folderList_tile)*e_val^2;
    d = 1+e/f;
    sample_size = round(a/d);
    neces_tiles = [neces_tiles; sample_size];
    tile_patch_per_2 = [tile_patch_per_2; 100*length(sub_folderList)/sample_size];
end
file_ID = folderName_clean_patch';
% contextual
varNames_contx_inter = {'Number of patches','Number of tiles', 'tiles excluded','needed patches', 'patch to tile ratio (%)', 'neces. patch to tile ratio (%)'};
cot_cohort_feat = [patches_cohort tiles_cohort exclude_tiles neces_tiles tile_patch_per tile_patch_per_2];
contextual_features = array2table(cot_cohort_feat,"VariableNames",varNames_contx_inter);
contextual_features = addvars(contextual_features,file_ID,'Before',1);
% Save the tables
writetable(contextual_features,[save_path_full cohort_name '_patches_information.csv']);
