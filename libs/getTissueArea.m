function a = getTissueArea( I )
%GETTISSUEAREA Summary of this function goes here
%   Detailed explanation goes here

J=rgb2gray(I);
%bw=im2bw(J,graythresh(J));
bw=~(J>220);
bw = bwareaopen(bw, 1000);
bw = bwareaopen(~bw, 100);
bw = ~bw;
%imshow(bw);
a=sum(sum(bw));

end

