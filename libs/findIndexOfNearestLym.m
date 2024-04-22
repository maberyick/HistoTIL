function [feat,featDiscription] = findIndexOfNearestLym(NucleiIndex,NucleiPoints,LymIndex,LymPoints)
Nuclei2LymDistV = [];
Nuclei2LymIndexV = [];%uniqueness of lym idx
Nuclei2LymDegreeV = [];
Nuclei2LymClosenessV = [];
Nuclei2LymBetweennessV = [];
Nuclei2LymEigenvectorV = [];
if isempty(LymPoints) || length(NucleiPoints)<1
    feat = zeros(1,7);
else
    for i = 1 : size(NucleiPoints,1)
        NucleiPoint =  NucleiPoints(i,:);
        N2LymDist = pdist2(NucleiPoint, LymPoints);
        [Dsorted,IdxLymByDist] = sort(N2LymDist);
        LymNearestPoint = LymPoints(IdxLymByDist(1),:);
        Nuclei2LymDistV = [Nuclei2LymDistV;Dsorted(1)];
        Nuclei2LymIndexV = [Nuclei2LymIndexV;IdxLymByDist(1)];
        
        if ~isempty(LymIndex.WeightedDegree)
            Nuclei2LymDegreeV = [Nuclei2LymDegreeV;LymIndex.WeightedDegree(IdxLymByDist(1))];
        end
        
        if ~isempty(LymIndex.WeightedCloseness)
            Nuclei2LymClosenessV = [Nuclei2LymClosenessV;LymIndex.WeightedCloseness(IdxLymByDist(1))];
        end
        
        if ~isempty(LymIndex.WeightedBetweenness)
            Nuclei2LymBetweennessV = [Nuclei2LymBetweennessV;LymIndex.WeightedBetweenness(IdxLymByDist(1))];
        end
        
        if ~isempty(LymIndex.WeightedEigenvector)
            Nuclei2LymEigenvectorV = [Nuclei2LymEigenvectorV;LymIndex.WeightedEigenvector(IdxLymByDist(1))];
        end
        
        %pdist2(NucleiPoint,LymNearestPoint)
        %plot([NucleiPoint(1),LymNearestPoint(1)],[NucleiPoint(2),LymNearestPoint(2)],'--')
        %drawArrow = @(x,y) quiver( x(1),y(1),x(2)-x(1),y(2)-y(1),0,'lineWidth',2,'MaxHeadSize',3);
        %     quiver(NucleiPoint(1),LymNearestPoint(1),NucleiPoint(2),LymNearestPoint(2),'MaxHeadSize',0.5)
        %drawArrow([LymNearestPoint(1),NucleiPoint(1)],[LymNearestPoint(2),NucleiPoint(2)]);
    end
    
    
    
    Nuclei2LymSumOfNormalizedCentralityIndex = sum([Nuclei2LymDegreeV./norm(Nuclei2LymDegreeV',1),Nuclei2LymClosenessV./norm(Nuclei2LymClosenessV',1),...
        Nuclei2LymBetweennessV./norm(Nuclei2LymBetweennessV',1),Nuclei2LymEigenvectorV./norm(Nuclei2LymEigenvectorV',1)],2);
    NuceliSelfSumofNormalizedCentralityIndex = sum([NucleiIndex.WeightedDegree./norm(NucleiIndex.WeightedDegree,1),NucleiIndex.WeightedCloseness./norm(NucleiIndex.WeightedCloseness,1)...,
        NucleiIndex.WeightedBetweenness./norm(NucleiIndex.WeightedBetweenness,1),NucleiIndex.WeightedEigenvector./norm(NucleiIndex.WeightedEigenvector,1)],2);
    
    % figure
    % plot(NucleiIndex.WeightedDegree/max(NucleiIndex.WeightedDegree),'*-'),hold on,title('WeightedDegree')
    % plot(Nuclei2LymDegreeV/max(Nuclei2LymDegreeV),'o-')
    %
    % figure
    % plot(NucleiIndex.WeightedCloseness/max(NucleiIndex.WeightedCloseness),'*-'),hold on,title('WeightedCloseness')
    % plot(Nuclei2LymClosenessV/max(Nuclei2LymClosenessV),'o-')
    %
    % figure
    % plot(NucleiIndex.WeightedBetweenness/max(NucleiIndex.WeightedBetweenness),'*-'),hold on,title('WeightedBetweenness')
    % plot(Nuclei2LymBetweennessV/max(Nuclei2LymBetweennessV),'o-')
    %
    % figure
    % plot(NucleiIndex.WeightedEigenvector/max(NucleiIndex.WeightedEigenvector),'*-'),hold on,title('WeightedEigenvector')
    % plot(Nuclei2LymEigenvectorV/max(Nuclei2LymEigenvectorV),'o-')
    %
    % figure
    % plot(NuceliSelfSumofNormalizedCentralityIndex/max(NuceliSelfSumofNormalizedCentralityIndex),'*-'),hold on,title('SumOfIndex')
    % plot(Nuclei2LymSumOfNormalizedCentralityIndex/max(Nuclei2LymSumOfNormalizedCentralityIndex),'o-')
    
    feat1 = norm(NucleiIndex.WeightedDegree/max(NucleiIndex.WeightedDegree)./Nuclei2LymDegreeV/max(Nuclei2LymDegreeV),1)/length(Nuclei2LymIndexV);
    %1-norm of the vector by ratio of nucleiDegree to LymDegree
    feat2 = norm(NucleiIndex.WeightedCloseness/max(NucleiIndex.WeightedCloseness)./Nuclei2LymClosenessV/max(Nuclei2LymClosenessV),1)/length(Nuclei2LymIndexV);
    %1-norm of the vector by ratio of nucleiCloseness to LymCloseness
    feat3 = norm(NucleiIndex.WeightedBetweenness/max(NucleiIndex.WeightedBetweenness)-Nuclei2LymBetweennessV/max(Nuclei2LymBetweennessV),1)/length(Nuclei2LymIndexV);
    %1-norm of the vector by nucleiBetweenness minus LymBetweenness
    feat4 = norm(NucleiIndex.WeightedEigenvector/max(NucleiIndex.WeightedEigenvector)./Nuclei2LymEigenvectorV/max(Nuclei2LymEigenvectorV),1)/length(Nuclei2LymIndexV);
    %1-norm of the vector by ratio of nucleiEigenvector to LymEigenvector
    
    feat5 = 0;
    if ~isempty(NuceliSelfSumofNormalizedCentralityIndex)
        feat5 = norm(NuceliSelfSumofNormalizedCentralityIndex/max(NuceliSelfSumofNormalizedCentralityIndex)./Nuclei2LymSumOfNormalizedCentralityIndex/max(Nuclei2LymSumOfNormalizedCentralityIndex),1)/length(Nuclei2LymIndexV);
    end
    %1-norm of the vector by ratio of nucleiSumIndex to LymSumIndex
    feat6 = length(unique(Nuclei2LymIndexV))/length(Nuclei2LymIndexV);
    %ratio of number of unique lym nodes to all lym nodes
    feat7 = sum(Nuclei2LymDistV)/length(Nuclei2LymIndexV);
    %average distance of nuclei to nearest lym
    
    % a = [NucleiIndex.WeightedDegree/max(NucleiIndex.WeightedDegree)-Nuclei2LymDegreeV/max(Nuclei2LymDegreeV)];
    % for i = 1:size(NucleiPoints,1)
    % drawArrow = @(x,y) quiver( x(:,1),y(:,1),x(:,2)-x(:,1),y(:,2)-y(:,1),0,'lineWidth',2,'MaxHeadSize',3);
    % q=drawArrow([LymPoints(Nuclei2LymIndexV(i),1),NucleiPoints(i,1)],[LymPoints(Nuclei2LymIndexV(i),2),NucleiPoints(i,2)]);
    % q.Color = [abs(a(i)) abs(a(i)) 1];
    % end
feat = [feat1,feat2,feat3,feat4,feat5,feat6,feat7];    
end

featDiscription = {'1-norm of the vector by ratio of nucleiDegree to LymDegree'...
    '1-norm of the vector by ratio of nucleiCloseness to LymCloseness'...
    '1-norm of the vector by nucleiBetweenness minus LymBetweenness'...
    '1-norm of the vector by ratio of nucleiEigenvector to LymEigenvector'...
    '1-norm of the vector by ratio of nucleiSumIndex to LymSumIndex'...
    'ratio of number of unique lym nodes to all lym nodes'...
    'average distance of nuclei to nearest lym'};