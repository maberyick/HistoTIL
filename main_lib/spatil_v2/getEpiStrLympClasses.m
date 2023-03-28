function classes = getEpiStrLympClasses(isLymphocyte,epiNuc)
%GETEPISTRLYMPCLASSES Summary of this function goes here
%   Detailed explanation goes here

classes=zeros(1,length(isLymphocyte));
classes(isLymphocyte & ~epiNuc)=1;
classes(isLymphocyte & epiNuc)=2;
classes(~isLymphocyte & ~epiNuc)=3;

end

