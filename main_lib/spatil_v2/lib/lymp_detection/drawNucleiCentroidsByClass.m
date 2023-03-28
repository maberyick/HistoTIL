function drawNucleiCentroidsByClass( image,centroids,class )

imshow(image);
hold on;

n=size(centroids,1);

for i=1:n
    if class(i)==1
        dot='g*';
    else
        dot='r*';
    end
    
    plot(centroids(i,1),centroids(i,2),dot);
end

hold off;
end