%       Pseudo-label Generation
%
%       This code is used to generate skeleton pseudo-labels for training ASMs.
%       Based on the probability map of the auroral regions output by the AAM,
%       the ridges of each luminous region are detected.
%
%       Version 1.0
%       Author: Qian Wang, Wangying Bai, Wei Zhang, Jinming Shi
%       Cited as: Automatically sketching Auroral Skeleton Structure in All-sky Image for Measuring Aurora Arcs
%
%
%
%
% %   Example:
%       % Read image from a local file
%       imdata = imread('ngc6543a.jpg');
%
%       % Read images from an ASI image file
%      imgFolder='F:\auroralData\labeled2003_38044\'; %File ASI image

%
%       % Read images from a file of probability map of aurora';
%       proFolder='.\segImgs\';% File: probability map of aurora';


% close all
clear all

SE2 = strel('disk', 1);

imgFolder='.\ASIImgs\'; %"F:\auroralData\labeled2003_38044\'; %File ASI image
proFolder='.\segImgs\';% File: probability map of aurora';

namelist=proFolder;
imgs=dir([namelist '*.bmp']);


iptsetpref('ImshowBorder','tight')

mask_in=zeros(440,440);
for i=1:440
    for j=1:440
        if (i-220)^2+(j-220)^2<180^2 %170 for 60 degree
            mask_in(i,j)=1;
        end
    end
end

polarT=zeros(440,440);
for i=1:440
    for j=1:440
        if (i-220)^2+(j-220)^2<151^2
            polarT(i,j)=12;
        elseif (i-220)^2+(j-220)^2<181^2
            polarT(i,j)=16;
        else
            polarT(i,j)=20;
        end
        
        
    end
end

SE2 = strel('diamond',2);
SE5 = strel('diamond',5);
H1=ones(7,7);
H2=ones(3,3);
H=fspecial('average',3);%fspecial('gaussian',hsize, sigma)??????

patchSize=71;
patchSizeHalf=floor(patchSize/2);
start_win=ceil(patchSize/2);
end_win=440-floor(patchSize/2);

HPatch_S=ones(31,31);
HPatch_L=ones(71,71);
imgSize=440;
lamda=1;
medfiltSize=5;

for i=1:size(imgs,1)
    
    %% read ASI image and probabiilty map of auroral area
    imgName=imgs(i).name;
    
    
    hh=str2num(imgName(:,11:12));
    mn=str2num(imgName(:,13:14));
    ss=str2num(imgName(:,15:16));
    UT=hh+mn/60.0+ss/3600.0;
    clear hh mn ss Timen
    
    
    img=double(imread([imgFolder imgName]));
    
    pMap=double(imread([proFolder imgName]));
    
    map=pMap;
    
    %% distance map
    map_fMed=medfilt2(map,[medfiltSize medfiltSize]);
    map_fMean=imfilter(map_fMed,H);
    
    if isempty(find(map_fMean>=100))
        continue
    end
    
    map_fMean(map_fMean<30)=0;
    bwMap = imbinarize(map_fMean, 'adaptive');
    bwMap2 = bwareafilt( bwMap,[100 inf]);% for connected area
    if isempty(find(bwMap2))
        continue
    end
    
    bwMap2=imfill(bwMap2,'holes');
    bwMap2=bwmorph(bwMap2,'spur',inf);
    
    edgeMap=edge(bwMap);
    [ASIcol,ASIrow]=meshgrid(1:440);
    ASIgrid=[ASIcol(:) ASIrow(:)];
    [edgerow,edgecol]=find(edgeMap);
    edgegrid=[edgecol edgerow];
    pd=pdist2(edgegrid,ASIgrid,'euclidean','Smallest',1);
    
    Lmiddle2=reshape(pd(1,:),440,440).* bwMap2;
    Lmiddle2= Lmiddle2./max( Lmiddle2(:));
    
    Lmiddle3=Lmiddle2.^lamda.* img.*pMap; %alpha parameter
    
    
    %% each area
    areas=regionprops(bwMap2,'Area','PixelIdxList','Image','Orientation');
    
    allRidge=zeros(440,440);
    figure
    imshow(img,[]);hold on
    for j=1:length(areas)
        if areas(j).Area>200
            layerMask=zeros(size(img));
            layerMask(areas(j).PixelIdxList)=1;
            DisMap_Layer=layerMask.*Lmiddle3;
            
            
            
            %% Option 1
            
%             DisMap_Layer_R=imrotate(DisMap_Layer,-areas(j).Orientation,'bilinear');
%             [fridge,~,lridge] = tfridge(DisMap_Layer_R,[1:length(DisMap_Layer_R)],0.01,'NumRidges',1,'NumFrequencyBins',10);
%             
%             P=zeros(size(DisMap_Layer_R));
%             for jj=1:length(fridge)
%                 P(fridge(jj)-1:fridge(jj)+1,jj)=255;
%             end
%             
%             ridge_Layer=imrotate(P,areas(j).Orientation);
%             ridge_Layer=bwmorph(ridge_Layer,'bridge');
%             [r,l]=size(ridge_Layer);
%             ridge_Layer_c= imcrop( ridge_Layer,[r/2-220  l/2-220  439 439]).*layerMask;
            
            %% Option 2
                        RidgeIn=ridge_AreaOnASI(UT,DisMap_Layer);
                        [r,l]=size(RidgeIn);
                        ridge_Layer_c= imcrop(RidgeIn,[r/2-220  l/2-220  439 439]).*layerMask;
                        [rr,cc]=find(ridge_Layer_c);
                        v=spcrv([cc,rr]',3);
                        plot(v(1,:),v(2,:),'r','Linewidth',2);hold on
            
            %%
                      
            
            allRidge= ridge_Layer_c+allRidge;
        end
       
        
        
        
    end
    %%  for show option 1.
%             figure
%             B = imoverlay(img./max(img(:)),allRidge,'red');
%             figure('visible','on')
%             imshow(B,[])
    
    
end
