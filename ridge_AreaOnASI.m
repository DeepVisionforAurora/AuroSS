function ridgeInArea=ridge_AreaOnASI(UT, GMap0)
%       Function: Detect ridge in an area
%      
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

angleNS=28.8664;
load ASImap2aacgm_G150.mat %²ÎÊý°üº¬ÓÐ£ºcgmLatG cgmLonG deltaMLT_G altitude;
load IntensityCorr4ASI_G150.mat %²ÎÊý°üº¬ÓÐ£ºIntensityCorrG;
cgmLat=cgmLatG;
cgmLon=cgmLonG;
deltaMLT=deltaMLT_G;
IntensityCorr=IntensityCorrG;
clear cgmLatG cgmLonG deltaMLT_G altitude IntensityCorrG


arcdate=cell(1,1);
aretime=cell(1,1);
nt=0;
ntpick=1;





m=size(GMap0,1);
n=size(GMap0,2);
r=220; %60 degree=164; 90:246
m1=-m/2:m/2-1;
n1=-n/2:n/2-1;
[x,y]=meshgrid(m1,n1);
circle=x.^2+y.^2;

circ_mask=zeros(m,n);
circ_mask(find(circle<=r*r))=1;
circ_mask(find(circle>r*r))=0;
GMap0=GMap0.*circ_mask;

GMap0=imrotate(GMap0,-61.1,'bilinear','crop');
GMap=zeros(512,512);
dx=x0-220;
dy=y0-220;
GMap(1+dy:440+dy,1+dx:440+dx)=GMap0;







%%--------------raw data to image---------
x=cgmLon; % longitude at each pixel
y=cgmLat; % latitude at each pixel
MLT=(UT+3)*ones(512,512)+deltaMLT;
%     [X,Y]=pol2cart(MLT*(2*pi/24)-pi/2,90-y);
% GMap=(Image1-bg_noise).*K;   % back noise and rescaling
% GMap=GMap.*IntensityCorr;
GMap1=GMap';
% Im=max(max(GMap(150:370,150:370)));
% c_Limn=[0 3000];
% I=mat2gray(GMap,c_Limn);
% clear c_Limn

%%---------------------------------------------------
%%---------------------------------------------------

[imagens0,imagex0,imagey0]=hatching1(GMap,angleNS-90,x0,y0); % vertical line of MN-MS
[mns,nns]=size(imagens0);
theta0=(1:nns)-ceil(nns/2);
theta0=90.0*theta0/theta0(1);% -90~90 equally divided
nline=161;         %121:´-60~60
r00=ones(1,nline)*r0;
for i=1:nline
    %         temp=abs(theta0-15*(i-4));
    temp=abs(theta0-1*(i-ceil(nline/2))); % theta0: degree at each pixel
    mintemp=min(temp);
    ind=min(find(temp==mintemp));% where is -60~60 on the MS-MN
    x00(i)=imagex0(ind); % where is -60~60 on the image
    y00(i)=imagey0(ind);
    r00(i)=r0*sqrt(90.0^2-(i-ceil(nline/2))^2)/90.0;
    clear temp mintemp ind
end
clear nns mns
imagens=zeros(nline,500);
imagex=zeros(nline,500);
imagey=zeros(nline,500);
indarc=zeros(nline,10);
arc=zeros(nline,10);
arctheta=zeros(nline,10);
ximage=zeros(nline,10);
yimage=zeros(nline,10);
xlon=zeros(nline,10);
ylat=zeros(nline,10);
Xlat=zeros(nline,10);
YMLT=zeros(nline,10);

for i=1:nline
    [a0,b0,c0]=hatching2(GMap,r00(i),x00(i),y00(i));% a0 is the intensity on MN-MS b0:x   c0:y
    [mns,nns(i)]=size(a0);
    imagens(i,1:nns(i))=a0;
    imagex(i,1:nns(i))=b0;
    imagey(i,1:nns(i))=c0;
    theta1=(1:nns(i))-ceil(nns(i)/2);
    theta1=90.0*theta1/theta1(1);
    magimage=a0;
    magimage=magimage-min(magimage(30:nns(i)-30));
    mmimage=max(magimage(10:nns(i)-10));
    magimage=1.0*magimage/mmimage;
    [a1,b1,c1,np(i)]=arcfinding(magimage,theta1);
    if np(i)>0   % number of peak
        indarc(i,1:np(i))=a1; %index of peak
        arc(i,1:np(i))=b1; % normalized value of peak
        arcn(i,1:np(i))=imagens(i,a1); % value of peak
        arctheta(i,1:np(i))=c1;
        for j=1:np(i)
            ximage(i,j)=imagex(i,a1(j)); % x of peak is signed to ith pararal line;
            yimage(i,j)=imagey(i,a1(j));
            xlon(i,j)=x(imagex(i,a1(j)),imagey(i,a1(j))); % longitude jing du
            ylat(i,j)=y(imagex(i,a1(j)),imagey(i,a1(j))); % latitude weidu
            YMLT(i,j)=MLT(imagex(i,a1(j)),imagey(i,a1(j)));
        end
    end
    if i==ceil(nline/2)
        xxx=theta1;
        yyy=a0;%magimage; red dot:
        arct0=c1;
    end
    %         xxx(i,1:nns(i))=theta1;
    %         yyy(i,1:nns(i))=magimage;
    clear mns mmimage magimage theta1 a0 b0 c0 a1 b1 c1
end


ntemp=ceil(nline/2);
if sum(np)~=0 % there must be a arc on MS-MN
    xpick=zeros(nline,np(ntemp));
    ypick=zeros(nline,np(ntemp));
    mltpick=zeros(nline,np(ntemp));
    

    
    for i=1:max(np)
        tempx0=ximage(ntemp,i);
        tempy0=yimage(ntemp,i);
        kpick=1;
        xpick(kpick,i)=xlon(ntemp,i);
        ypick(kpick,i)=ylat(ntemp,i);
        mltpick(kpick,i)=YMLT(ntemp,i);
        nwarn=0;
        for j=1:ntemp-1 % forward
            if np(ntemp-j)>0 & nwarn<nline/2
                tempx1=ximage(ntemp-j,:); % x of peaks
                tempy1=yimage(ntemp-j,:); % y of peaks
                dtemp=abs(tempx1-tempx0)+abs(tempy1-tempy0);
                mindtemp=min(dtemp);
                if mindtemp<100 % find the point in a line
                    ind=1;%min(find(dtemp==mindtemp));
                    tempx0=tempx1(ind); % move to next point on line
                    tempy0=tempy1(ind);
                    kpick=kpick+1;
                    xt001(kpick,i)=ximage(ntemp-j,ind);
                    yt001(kpick,i)=yimage(ntemp-j,ind);
                    xpick(kpick,i)=xlon(ntemp-j,ind); % jingdu of next point
                    ypick(kpick,i)=ylat(ntemp-j,ind);
                    mltpick(kpick,i)=YMLT(ntemp-j,ind);
                    clear ind
                else %not on a line
                    nwarn=nwarn+1;
                end
                clear tempx1 tempy1 dtemp mindtemp
            else
                nwarn=nwarn+1;
            end
        end
        
        tempx0=ximage(ntemp,i);
        tempy0=yimage(ntemp,i);
        nwarn=0;
        for j=1:ntemp-1
            if np(ntemp+j)>0 & nwarn<nline/2
                tempx1=ximage(ntemp+j,:);
                tempy1=yimage(ntemp+j,:);
                dtemp=abs(tempx1-tempx0)+abs(tempy1-tempy0);
                mindtemp=min(dtemp);
                if mindtemp<100
                    ind=1;
%                     ind=min(find(dtemp==mindtemp));
                    tempx0=tempx1(ind);
                    tempy0=tempy1(ind);
                    kpick=kpick+1;
                    xt001(kpick,i)=ximage(ntemp+j,ind);
                    yt001(kpick,i)=yimage(ntemp+j,ind);
                    xpick(kpick,i)=xlon(ntemp+j,ind);
                    ypick(kpick,i)=ylat(ntemp+j,ind);
                    mltpick(kpick,i)=YMLT(ntemp+j,ind);
                    clear ind
                else
                    nwarn=nwarn+1;
                end
                clear tempx1 tempy1 dtemp mindtemp
            else
                nwarn=nwarn+1;
            end
        end
        clear tempx0 tempy0 kpick nwarn
    end
    
    

        
        ridgeInArea=zeros(512,512);
        xt001_1=xt001(xt001~=0);
        yt001_1=yt001(yt001~=0);
        
        for i=1:length( xt001_1)
            ridgeInArea( yt001_1(i), xt001_1(i))=255;
        end
        ridgeInArea=imrotate(ridgeInArea,61.1,'crop');
    
%         figure;
%         imshow(GMap,[]);hold on
%          for i=1:length( xt001_1)
%             plot(  xt001_1(i),yt001_1(i),'ro'); hold on
%         end
        
%     end
end

