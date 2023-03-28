function [features] = getSpaTILFeaturesv2(lympCentroids_cur,nonLympCentroids_cur)
%GETSPATILFEATURES V2 Gets a set of spatial quantitative features (SpaTIL) V2
% relating to the spatial architecture of TILs, co-localization of TILs and
% cancer nuclei, and density variation of TIL clusters. 
coords={lympCentroids_cur,nonLympCentroids_cur,};
[spaFeat_2gr_37,~]=getSpaTILFeatures_v2(coords);
[spaFeat_2gr_40,~]=getSpaTILFeatures_v2(coords);
[spaFeat_2gr_42,~]=getSpaTILFeatures_v2(coords);
[spaFeat_2gr_43,~]=getSpaTILFeatures_v2(coords);
% combine features
features = [spaFeat_2gr_37 spaFeat_2gr_40 spaFeat_2gr_42 spaFeat_2gr_43];
end