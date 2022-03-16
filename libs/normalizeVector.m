function arr=normalizeVector(vect,inverted)
%NORMALIZEVECTOR maps a vector to the range [0,1]
if inverted==1 % major = black
    arr = ((max(vect) - vect) / (max(vect) - min(vect)));
else
    arr = ((vect - min(vect)) / (max(vect) - min(vect)));
end
 

end