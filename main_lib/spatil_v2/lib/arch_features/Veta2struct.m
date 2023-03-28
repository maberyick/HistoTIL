function [output] = Veta2struct(nuclei, centroids)
% nuclei=nuclei;

%Nuclei come in as a cell
%Centroids come in as an array

%Output is a cell
if ~isempty(cellfun(@isempty,nuclei))
    if max(cellfun(@isempty,nuclei))==0
        for i=1:1:length(nuclei)
            clear subsetnuclei
            subsetnuclei=cell2mat(nuclei(i));
            output(i).r=subsetnuclei(:,1)';
            output(i).c=subsetnuclei(:,2)';
            output(i).centroid_r=centroids(i,2);
            output(i).centroid_c=centroids(i,1);
        end
    else
        for i=1:1:length(nuclei)
            output(i).r=[];
            output(i).c=[];
            output(i).centroid_r=[];
            output(i).centroid_c=[];
        end
    end
else
            output(1).r=[];
            output(1).c=[];
            output(1).centroid_r=[];
            output(1).centroid_c=[];
end
    