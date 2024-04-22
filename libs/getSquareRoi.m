function R = getSquareRoi( I,centroid,roiW,roiH )
%GETSQUAREROI extracts a square region of size roiW x roiH around the
%centroid.
[lim,~]=size(I);
len=max(roiW,roiH)/2;
x1=round(centroid(1)-len);
x2=round(centroid(1)+len);
y1=round(centroid(2)-len);
y2=round(centroid(2)+len);

w=x2-x1;
h=y2-y1;

if w>h
    y2=y2+1;
elseif h<w
    x2=x2+1;
end

R=uint8(zeros(x2-x1+1,x2-x1+1));

if x1<1
    x1=1;
end
if y1<1
    y1=1;
end

if y2>lim
    y2=lim;
end

if x2>lim
    x2=lim;
end

R(1:y2-y1+1,1:x2-x1+1)=I(y1:y2,x1:x2);

end
