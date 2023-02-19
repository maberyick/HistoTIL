function [full_feat_set] = stats_calc_feature(file_comp)
% prepare the variables and cases
%file_list = {};
%% calculate the basic descriptive statistic
% https://www.mathworks.com/help/matlab/data_analysis/descriptive-statistics.html
% max: Maximum value
value_bds_1 = max(file_comp);
% mean: Average or mean value
value_bds_2 = mean(file_comp);
% median: Median value
value_bds_3 = median(file_comp);
% min: Smallest value
value_bds_4 = min(file_comp);
% mode: Most frequent value
value_bds_5 = mode(file_comp);
% std: Standard deviation
value_bds_6 = std(file_comp);
% var: Variance, which measures the spread or dispersion of the values
value_bds_7 = var(file_comp);
% More basics
% https://www.mathworks.com/help/matlab/descriptive-statistics.html
% rms: Root-mean-square value (Basic)
value_bds_8 = rms(file_comp);
% iqr: Interquartile range of data set (Percentile and Quantile)
value_bds_9 = iqr(file_comp);
%% calculate the advanced descriptive statistic: Central Tendency and Dispersion
% https://www.mathworks.com/help/stats/descriptive-statistics.html
% geomean: Geometric mean
try
    value_ctd_1 = geomean(file_comp);
catch
    %disp('geoamtric mean: it did not contain negative or complex values. return 0s')
    [~,b] = size(file_comp);
    value_ctd_1 = zeros(1,b);
end
% harmmean: Harmonic mean
value_ctd_2 = harmmean(file_comp);
% trimmean:	Mean, excluding outliers at 10% of the input values
value_ctd_3 = trimmean(file_comp, 10);
% kurtosis:	Kurtosis
value_ctd_4 = kurtosis(file_comp);
% moment: Central moment with the order of 2 (variance)
value_ctd_5 = moment(file_comp,2);
% skewness:	Skewness
value_ctd_6 = skewness(file_comp);
%
%% calculate the advanced descriptive statistic: Range, Deviation, and Z-Score
% https://www.mathworks.com/help/stats/descriptive-statistics.html
% range: Range of values
value_rdz_1 = range(file_comp);
% mad: Mean or median absolute deviation
value_rdz_2 = mad(file_comp);
% zscore:	Standarized z-score and the mean value calculated from it
value_rdz_3 = mean(zscore(file_comp));
%
%% Put the feature set together
full_feat_set = [value_bds_1 value_bds_2 value_bds_3 value_bds_4 ...
    value_bds_5 value_bds_6 value_bds_7 value_bds_8 value_bds_9 ...
    value_ctd_1 value_ctd_2 value_ctd_3 value_ctd_4 value_ctd_5 ...
    value_ctd_6 value_rdz_1 value_rdz_2 value_rdz_3];