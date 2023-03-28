function [bounds] = DLMask2bounds_function(originalpath, maskpath)
% debug=1;
% originalpath='F:\Jon\Test images\Nuclei ground truth\10256_500_f00001_original.tif';
% groundtruthpath='F:\Jon\Test images\Nuclei ground truth\10256_500_f00001_mask.png';
% originalpath='F:\Jon\Breast Cancer Data\DL Results\Images';
% files = dir(originalpath);
% j=1;
% k=1;

%This expects all mask files to have the word mask in the title, and only
%original pictures to be .tifs
% for i=1:length(files)
%     if ~isempty(strfind(files(i).name, 'mask'));
%         filepath(j).mask=strcat(originalpath, '\', files(i).name);
%         j=j+1;
%     elseif ~isempty(strfind(files(i).name, 'tif'));
%         filepath(k).original=strcat(originalpath, '\', files(i).name);
%         k=k+1;
%     end
% end
if ischar(originalpath)
    img=imread(originalpath);
else
    img=originalpath;
end
if ischar(maskpath)
    mask=imread(maskpath);
else
    mask=maskpath;
end


[originallength, originalwidth, ~]=size(img);
mask=imresize(mask, [originallength, originalwidth]);

if ~islogical(mask)
    if max(max(mask)) > 1
        nucleibw=im2bw(mask, .4);
    else
        nucleibw=logical(mask);
    end
else
    nucleibw=mask;
end
nucleibw=imfill(nucleibw, 'holes');
% nucleibw = bwareaopen(nucleibw, 50);
% figure;
% imshow(nucleibw);
% title('before small object removal, 0.4 threshold');

bounds=bwnuclei2bounds(nucleibw);





% figure;
% imshow(img);
% axis image;
% %     axis off;    
% hold on;
% [lengthDL, widthDL]=size(bounds);
% 
%     for k = 1:widthDL
%         plot(bounds(k).c, bounds(k).r, 'g-', 'LineWidth', 2);
%     end
% title('Deep Learning segmentation');
end