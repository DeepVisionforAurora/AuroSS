%        Calculate the arc tilt to measure auroral orientation
%        AuroSS in ASI image are mapped to MLT-MLAT coordinates
%
%       This code is used to calculate the arc tilt to measure auroral orientation.    
%       The obtained AuroSS  is converted into geomagnetic coordinate reference frame.  
%       We calculate the arc tilt to measure auroral orientation.
%   
%       Version 1.0
%       Author: Qian Wang, Wangying Bai, Wei Zhang, Jinming Shi
%       Cited as: Automatically sketching Auroral Skeleton Structure in All-sky Image for Measuring Aurora Arcs
% 
% 
% 
% 
% %   Example:
%       % Name of AuroSS map
%       filename='N20031221G081212.bmp'; % AuroSS name
%       Data naming convention:
%       Each auroral data has a unique file name that contains the N/S, date, band, and time information.
%       Filename example: N20031221G081212.
%           N: North pole
%           20031221: December 21, 2003
%       	G: G-band
% 	        081212: 8:12:12 UT”

%
%       % Read AuroSS map from a local file
%       SMap0=double(imread(['.\skeImgs\'  filename]));




clear;clc;
close all
Band='G';

band='557.7nm';
bg_noise=ones(512,512)*564;
c_Lim=[0 5000];
%             c_LimRayleigh=[0 1500];
c_LimRayleigh=[0 3000];
K=1.0909;

x0=261;
y0=257;
r0=246; % radius
h=150;
angleNS=28.8664;
load ASImap2aacgm_G150.mat %cgmLatG cgmLonG deltaMLT_G altitude;
load IntensityCorr4ASI_G150.mat %IntensityCorrG;
cgmLat=cgmLatG;
cgmLon=cgmLonG;
deltaMLT=deltaMLT_G;
IntensityCorr=IntensityCorrG;
clear cgmLatG cgmLonG deltaMLT_G altitude IntensityCorrG


arcdate=cell(1,1);
aretime=cell(1,1);
nt=0;


filename='N20040116G050513.bmp'; % AuroSS name
% [Image1,Date1,Time1,Tag1,Exposure1]=OpenImg2004(filename);
SMap0=double(imread(['.\skeImgs\'  filename]));
% a=load(['.\psudoL_bad\disMap-N20031222G042701.mat']);
% SMap0=a.Lmiddle2;

m=440;
n=440;
r=170; %60 degree=164; 90:246
m1=-m/2:m/2-1;
n1=-n/2:n/2-1;
[x,y]=meshgrid(m1,n1);
circle=x.^2+y.^2;

circ_mask=zeros(m,n);
circ_mask(find(circle<=r*r))=1;
circ_mask(find(circle>r*r))=0;
SMap0=SMap0.*circ_mask;





SMap0=imrotate(SMap0,-61.1,'bilinear','crop');
SMap=zeros(512,512);
dx=x0-220;
dy=y0-220;
SMap(1+dy:440+dy,1+dx:440+dx)=SMap0;





% Timen=char(Time1);
hh=str2num(filename(:,11:12));
mn=str2num(filename(:,13:14));
ss=str2num(filename(:,15:16));
UT=hh+mn/60.0+ss/3600.0;
clear hh mn ss Timen

%%---------------------
x=cgmLon; % longitude at each pixel
y=cgmLat; % latitude at each pixel
MLT=(UT+3)*ones(512,512)+deltaMLT;
%     [X,Y]=pol2cart(MLT*(2*pi/24)-pi/2,90-y);


ntpick=1;
% GMap=(Image1-bg_noise).*K;   % È¥³ýCCD±³¾°ÔëÒô,¼ÆËã¼«¹âÇ¿¶È
% GMap=GMap.*IntensityCorr;
SMap1=SMap';
% Im=max(max(GMap(150:370,150:370)));
% c_Limn=[0 3000];
% I=mat2gray(GMap,c_Limn);
% clear c_Limn


%%  find intersection
[a0,b0,c0]=hatching2(SMap,r0,x0,y0);%
npink_MNMS=size(a0,2);
image_MNMS=a0;
imagex_MNMS=b0;
imagey_MNMS=c0;
theta1_s=(1:npink_MNMS)-ceil(npink_MNMS/2);
theta1_S=90.0*theta1_s/theta1_s(1); %-90:90

image_MNMS=image_MNMS-min(image_MNMS(30:end-30));
mmimage=max(image_MNMS(10:end-10));
image_MNMS=image_MNMS/mmimage;
[pks,locs] = findpeaks(image_MNMS);

if length(pks)>0
    ximage_s=imagex_MNMS(locs);%ximage: peak; imagex: all on line
    yimage_s=imagey_MNMS(locs);
    for ipick_MNMS=1:length(locs)
        xlon_MNMS(ipick_MNMS)=x(imagex_MNMS(locs(ipick_MNMS)),imagey_MNMS(locs(ipick_MNMS)));
        ylat_MNMS(ipick_MNMS)=y(imagex_MNMS(locs(ipick_MNMS)),imagey_MNMS(locs(ipick_MNMS)));
        YMLT_MNMS(ipick_MNMS)=MLT(imagex_MNMS(locs(ipick_MNMS)),imagey_MNMS(locs(ipick_MNMS)));
    end
end

[Xintersection_MNMS,Yintersection_MNMS]=pol2cart(YMLT_MNMS*(2*pi/24)-pi/2,90-ylat_MNMS); % weidu: ypick


XYpick_line(:,1)=Xintersection_MNMS';
XYpick_line(:,2)=Yintersection_MNMS';



%% ---------------------------------------------------
bw_img1 = im2bw(SMap, 0.8); %binarization
bw_img1 = bwareafilt(bw_img1,[200 inf]);

bw2=imfill(bw_img1,'holes');
img_seg=bwmorph(bw2,'spur',inf);
out = bwskel(img_seg,'MinBranchLength',20);
regionsInSke=regionprops(out,'Area','PixelList');
clear skeLength skeOrientation





numArc=length(regionsInSke);

skeLength=cat(1,regionsInSke.Area);
maxSkeLength=max(skeLength);

ximage_s=zeros(maxSkeLength,numArc);
yimage_s=zeros(maxSkeLength,numArc);
xlon_s=zeros(maxSkeLength,numArc);
ylat_s=zeros(maxSkeLength,numArc);
MLT_s=zeros(maxSkeLength,numArc);
for numSke=1:numArc
    
    
    ximage_s(1:skeLength(numSke),numSke)=regionsInSke(numSke).PixelList(:,1);
    yimage_s(1:skeLength(numSke),numSke)=regionsInSke(numSke).PixelList(:,2);
    
    
    
    
    for j=1:skeLength(numSke)
        xlon_s(j,numSke)=x(regionsInSke(numSke).PixelList(j,1),regionsInSke(numSke).PixelList(j,2));
        ylat_s(j,numSke)=y(regionsInSke(numSke).PixelList(j,1),regionsInSke(numSke).PixelList(j,2));
        MLT_s(j,numSke)=MLT(regionsInSke(numSke).PixelList(j,1),regionsInSke(numSke).PixelList(j,2));
    end
    ind=find(xlon_s(:,numSke)~=0);
    parc_s=polyfit(xlon_s(ind,numSke),ylat_s(ind,numSke),1);
    xarc_s=xlon_s(1:skeLength(numSke),numSke);
    yarc_s=polyval(parc_s,xarc_s);
    
    
    
    
    % %     tx = txsite('Name','MathWorks','Latitude',yarc_s(1),'Longitude',xarc_s(1));
    % %     rx = rxsite('Name','Fenway Park','Latitude',yarc_s(end),'Longitude',xarc_s(end));
    % %     angle(tx,rx)
    % %
    
    
    mltarc_s=(xarc_s-x(255,255))/15.0+MLT(255,255);
    [Xarc_s,Yarc_s]=pol2cart(mltarc_s*(2*pi/24)-pi/2,90-yarc_s);
    Xarcn(numSke)={Xarc_s}; % for plot
    Yarcn(numSke)={Yarc_s};
    
    
    
    XYarc_line(:,1)=Xarc_s;
    XYarc_line(:,2)=Yarc_s;
    pd=pdist2(XYarc_line,XYpick_line,'euclidean','Smallest',1);
    [v,index]=min(pd);
    Xsec(numSke)=Xintersection_MNMS(index);
    Ysec(numSke)=Yintersection_MNMS(index);
    
    yptt=spline(Xarc_s,Yarc_s);
    dpyy=fnder(yptt);
    slopep(numSke)=ppval(dpyy,Xsec(numSke)); %%%%%%%%%%%%%%%%%%%%%%%%
    arc_orientation(numSke)=atand(slopep(numSke));
    
       
    tilttemp=(YMLT_MNMS(numSke)-6)*15-arc_orientation(numSke)-90;
    if tilttemp<-90
        tilttemp=tilttemp+180;
    end
    if tilttemp>90
        tilttemp=tilttemp-180;
    end
    arctilt(numSke)=tilttemp;
    arcmlt(numSke)=YMLT_MNMS(numSke);
    
    
    clear ytt dyy  XYarc_line
end


[Xpick_s,Ypick_s]=pol2cart(MLT_s*(2*pi/24)-pi/2,90-ylat_s);



%%  plot ---------------------------------------------------
figure
imshow(SMap,[])

figure
hndlel=pcolor(x,y,SMap1); % Intensity Correct
shading interp;
colormap(gray);
%         set(gca,'CLim',c_LimRayleigh);   %   É«±ê³ß¶È
hold on

%     plot(xpick,ypick,'og','markersize',2)
%     plot(xlon_s,ylat_s,'.g','markersize',2)
%

%     plot(xarc,yarc,'r','linewidth',1)
xlabel('Mag Lon (deg)','fontsize',14,'fontname','Times New Roman')
ylabel('Mag Lat (deg)','fontsize',14,'fontname','Times New Roman')
set(gca,'xlim',([90,140]),'xTick',[100:10:130])
set(gca,'ylim',([70,82]),'yTick',[72:2:80])
set(gca,'xminortick','on','yminortick','on','tickdir','out')
set(gca,'Box','on')
set(gca,'fontsize',14,'fontname','Times New Roman')
%         axis off
%         text(112,81,'N','color','g','fontsize',18,'fontname','Times New Roman')
%         text(133,76,'E','color','g','fontsize',18,'fontname','Times New Roman')
printname='arcalignt13';
set(gcf,'PaperType','A4');
print(gcf,'-dtiff',printname);
%         print(gcf,'-dpng',printname);
%         close(gcf);


figure

%         ax3=axes('position',[0.4500 0.330 0.510 0.340]);
half_Axis4AACGM1;
hold on
%         [X,Y]=pol2cart(MLT*(2*pi/24)-pi/2,90-y);
%         hndlel=pcolor(X,Y,GMap1); % Intensity Correct
%         shading interp;
plot(Xpick_s,Ypick_s,'og','markersize',5)

plot(Xintersection_MNMS,Yintersection_MNMS,'.r','markersize',18)
yt=(-1:30)+0.5;
for i=1:size(ylat_s,2)
    plot(Xarcn{i},Yarcn{i},'r','linewidth',2)
    xt1=(yt-Ysec(i))/slopep(i)+Xsec(i);
    plot(xt1,yt,'-b','linewidth',1)
end

clear slopepn arc_orientation XYarc_line  arcmlt arctilt





clear GMap Gmap1
clear Image1 Date1 Time1 Tag1 Exposure1






