function [epcellout, stcellout, labels]=epstroma_features_2(img, epitheliummask, stromamask, lpath, i)
% lpath = trimfilepath;
datawidth=267;
m=0;
% OCMfeature=collage_all_funct(img(:,:,1), '');

% [~, datawidth]=size(lowdatabase);
ep_objects = bwconncomp(epitheliummask);
if ~isempty(ep_objects.PixelIdxList)
%         [epbounds]=get_rid_of_short_bounds(epbounds, 10);
        [epfeats,descriptionep] = extract_all_features_Jon(epitheliummask, img);
        %need remove edge nuclei
        t1=tic;
%         [epcolfeats, epcollabels]=collage_cut_fromtotal(OCMfeature, epitheliummask, '');
%         [epcolfeats, epcollabels]=collage_nuclei_funct(img(:,:,1), epitheliummask, '');
%         disp([num2str(toc(t1)/60), 'minutes']);
%         epfeatstotal=[epfeats{:},epcolfeats];
%             eplabels=[descriptionep, epcollabels];
        epfeatstotal=[epfeats{:}];
        eplabels=[descriptionep];
        m=1;%check to see if labels holds useful info
        [~, datawidth]=size(epfeatstotal);
else
    epfeatstotal=cell(1,datawidth);
end
st_objects = bwconncomp(stromamask);        
if ~isempty(st_objects.PixelIdxList)       
%         [stbounds]=get_rid_of_short_bounds(stbounds, 10);
        [stfeats,descriptionst] = extract_all_features_Jon(stromamask,img);
%         [stcolfeats, stcollabels]=collage_cut_fromtotal(OCMfeature, stromamask, '');
%         [stcolfeats, stcollabels]=collage_nuclei_funct(img(:,:,1), stromamask, '');
%         stfeatstotal=[stfeats{:},stcolfeats];
%         stlabels=[descriptionst, stcollabels];
        stfeatstotal=[stfeats{:}];
        stlabels=[descriptionst];
        [~, datawidth]=size(stfeatstotal);
        
        m=2; %check to see if labels holds useful info
else
    stfeatstotal=cell(1,datawidth);
end


datawidth=datawidth+3;
%Get the labels written correctly
if m==0
    labels=cell(1,datawidth);
elseif m==1
    labels=eplabels;         
elseif m==2
    labels=stlabels;
end


labels{1, datawidth-2}='Ep=0,stroma=1';
labels{1, datawidth-1}='Name';
labels{1, datawidth}='ODX Score';
epcellout=featurecell2array(epfeatstotal,labels);
epcellout{2,datawidth-2}=0;
epcellout{2,datawidth-1}=lpath(i).originalname;
epcellout{2, datawidth}=lpath(i).odx;

stcellout=featurecell2array(stfeatstotal,labels);
stcellout{2, datawidth-2}=1;
stcellout{2, datawidth-1}=lpath(i).originalname;
stcellout{2, datawidth}=lpath(i).odx;

% for k1=1:datawidth
%             lowdatabase{i+1,k1}=epcellout{k1};
%             lowdatabase{i+2,k1}=stcellout{k1};
%         end
%     end