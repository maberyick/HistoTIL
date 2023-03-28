function plotClusters(image,coords,M,args)

%% Processing args

if ~isfield(args,'printGroupId')
    args.printGroupId=false;
end

if ~isfield(args,'printCenter')
    args.printCenter=false;
end

if ~isfield(args,'fillAlpha')
    args.fillAlpha=.25;
end

if ~isfield(args,'lineWidth')
    args.lineWidth=4;
end

if ~isfield(args,'lineAlpha')
    args.lineAlpha=.5;
end

if ~isfield(args,'useBlackLine')
    args.useBlackLine=true;
end

if ~isfield(args,'plotPoints')
    args.plotPoints=false;
end

if ~isfield(args,'plotDelaunay')
    args.plotDelaunay=false;
end

if ~isfield(args,'plotVoronoi')
    args.plotVoronoi=false;
end

if ~isfield(args,'plotCenterGrouping')
    args.plotCenterGrouping=false;
end

printGroupId=args.printGroupId;
fillAlpha=args.fillAlpha;
lineWidth=args.lineWidth;
lineAlpha=args.lineAlpha;
useBlackLine=args.useBlackLine;
plotPoints=args.plotPoints;
printCenter=args.printCenter;
plotVoronoi=args.plotVoronoi;
plotDelaunay=args.plotDelaunay;
plotCenterGrouping=args.plotCenterGrouping;

%%

colors={'b','y','g','r','c','m','w'};

numGroups=length(coords);

imshow(image);
hold on;
for k = 1:numGroups
    
    matrix=M{k};
    nodes=coords{k};
    
    [~,~,clusters] = networkComponents(matrix);
    clusters(cellfun(@length,clusters)<=2)=[];
    
    numClusters=length(clusters);
    
    cent=[];
    
    for i = 1:numClusters
        member = clusters{i};
        
        
        col = nodes(member,1);
        row = nodes(member,2);
        [pol,a] = convhull(col,row);
        
        %if a>30000
            x=col(pol);
            y=row(pol);
            cent=[cent;[mean(x) mean(y)]];
            
            s=plot(x,y,colors{k},'LineWidth',lineWidth);
            s.Color(4)=lineAlpha;
            
            if fillAlpha>0
                f=fill(x,y,colors{k},'LineStyle','none');
                set(f,'facealpha',fillAlpha)
            end
            
            if useBlackLine
                plot(x,y,'k','LineWidth',1);
            end
        %end
    end
    
    if ~isempty(cent)
        if printGroupId
            text(cent(:,1),cent(:,2),num2str(k),'Color','black','FontSize',20);
        end
        
        if printCenter
            scatter(cent(:,1),cent(:,2),100,[colors{k} 'o'],'filled');
        end
        
        if plotPoints
            scatter(nodes(:,1),nodes(:,2),100,[colors{k+1} 'o'],'filled');
        end
        
        if plotDelaunay
            DT=delaunay(cent(:,1),cent(:,2));
            triplot(DT,cent(:,1),cent(:,2),'c-','LineWidth',3);
        end
        
        if plotVoronoi
            [vx,vy]=voronoi(cent(:,1),cent(:,2));
            plot(cent(:,1),cent(:,2),'g*',vx,vy,'g-','LineWidth',5);
        end
        
        
        if plotCenterGrouping
            
            vect=getSumNodeWeightsThreshold( cent,'euclidean',.0005);
            vect=normalizeVector(vect);
            
            for i=1:size(cent,1)
                    hsl=getColorHeatMap(vect(i));
                    rgb=hsl2rgb(hsl);
                    plot(cent(i,1),cent(i,2),'or',...
                        'MarkerSize',15,...
                        'MarkerEdgeColor',[0,0,0],...
                        'MarkerFaceColor',[rgb(1,1,1),rgb(1,1,2),rgb(1,1,3)]);%,...
                    %text(nucleiCentroids(i,1),nucleiCentroids(i,2), sprintf('%d',i), 'Color', 'b')
            end
        end
        
        
    end
    
    
    
    
    
end
hold off;

end