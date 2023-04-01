function feat_extraction_main(varargin)
% Paths
addpath(genpath('C:\Users\cbarr23\Documents\HistoTIL\'))
%addpath(genpath('/home/maberyick/CCIPD_Research/Github/HistoTIL'))
folder_matpatches = varargin{1};
folder_pyepistroma = varargin{2};
folder_matcellmask = varargin{3};
folder_savepath = varargin{4};
folder_type = varargin{5};

% tmp linux
%folder_matpatches = '/home/maberyick/CCIPD_Research/Github/HistoTIL/testimage/patch/';
%folder_pyepistroma = '/home/maberyick/CCIPD_Research/Github/HistoTIL/testimage/epimask/';
%folder_matcellmask = '/home/maberyick/CCIPD_Research/Github/HistoTIL/testimage/nucleimask/';
%folder_savepath = '/home/maberyick/CCIPD_Research/Github/HistoTIL/testimage/histotil_features/';
%folder_type = 'tile_type';

% tmp windows
folder_matpatches = 'D:\Projects\SCLC\roswell\batch_b\tiles\';
folder_pyepistroma = 'D:\Projects\SCLC\roswell\batch_b\epimask\';
folder_matcellmask = 'D:\Projects\SCLC\roswell\batch_b\nucleimask\';
folder_savepath = 'D:\Projects\SCLC\roswell\batch_b\histotil_features\';
folder_type = 'folder_type';

% Quality parameters for Nuclear and peri-nuclear features
tileSize = 2048; % Size of the patch generated by the deep learning models
quality.redChannelLim=100; % ignore the current tile if the mean intensity of R channel is less than this threshold
quality.tissueLim=.1; % process the current tile only if it contains at least this percentage of tissue
quality.blurLimit=.2; % process the current tile only if its blurriness score is below this value
quality.saturationLim=tileSize*20;
quality.imgarea = 0.2*(tileSize*tileSize);
feat_extraction_wsi(folder_matpatches,folder_pyepistroma,folder_matcellmask,folder_savepath,quality,folder_type)
end