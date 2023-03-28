function stats = getFeatureStats(featVector)
%GETFEATURESTATS Summary of this function goes here
%   Detailed explanation goes here

 [r,c]=size(featVector);
 if r==0
     stats = nan(1,c*8);
 else
    stats=[nansum(featVector,1),nanmean(featVector,1),nanstd(featVector,0,1),...
        nanmedian(featVector,1),max(featVector,[],1),min(featVector,[],1),...
        kurtosis(featVector,1,1),skewness(featVector,1,1) ];
end

end

