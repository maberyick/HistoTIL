function [centroids,polygons]=getClusterCentroids(members,nodes)
numMembers=length(members);
centroids=[];
polygons={};
counter=1;
for i = 1:numMembers
    member = members{i};
    col = nodes.centroid_c(member)';
    row = nodes.centroid_r(member)';
    numNodes=length(col);
    if  numNodes >  2
        [k,~] = convhull(col,row);
        cx=mean(col(k));
        cy=mean(row(k));
        centroids=[centroids;[cx cy]];
        polygons{counter}=[col(k) row(k)];
        counter=counter+1;
    end
end

end