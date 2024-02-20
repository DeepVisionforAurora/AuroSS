function [indarc,arc,arctheta,np]=arcfinding(image,theta)

simage=smooth(image,5);        % 先作个平滑，把一些小波动去掉
ind=find(abs(theta)<70);       % 不考虑边界附近的数据
image1=simage(ind);
[m,n]=size(image1);           % m为数组个数，n=1
peak0=zeros(10,1);
indpeak0=zeros(10,1);          % 认为弧的数量（峰数）不多于10个
np=0;                          % 计数，弧的个数
for i=2:(m-1)
    np0=0;                 % sign：0表示没有满足条件的峰存在
    % -------------峰判断： 比两边都大；峰谷值差大于0.1。-------------
    if (image1(i)>image1(i-1)) & (image1(i)>=image1(i+1))           
        peak1=image1(i);
        j=i-2;
        while (j>0) & (image1(j)<=peak1 & np0==0)          
            di=peak1-image1(j);
            j=j-1;
            if di>0.1                                       % 峰左侧存在比其小0.1的值
                np0=1;
                j=0;
            end
            clear di
        end
        if np0==1
            np0=0;
            j=i+2;
            while (j<=m) & (image1(j)<=peak1 & np0==0)
                di=peak1-image1(j);
                j=j+1;
                if di>0.1
                    np0=1;
                    j=m+1;
                end
                clear di
            end
        end
        clear j
    end
    if np0==1                   % 如果存在满足条件的峰，找到其位置和峰值，并让计数加1
        np=np+1;
        peak0(np)=peak1;
        indpeak0(np)=i+ind(1)-1;
    end
    clear np0
end
   

if np~=0
    indp=find(peak0>0);
    peak=peak0(indp);
    indpeak=indpeak0(indp);
    ptheta=theta(indpeak);    % 把非0的弧赋值到新数组
    clear indp peak0 indpeak0
    if np>1
        i=2;
        while i<=np
            valley=min(simage(indpeak(i-1):indpeak(i)));
            if peak(i)<=peak(i-1)
                if peak(i)-valley<peak(i)/4.0
                    if np>2
                        peak(i:(np-1))=peak((i+1):np);
                        indpeak(i:(np-1))=indpeak((i+1):np);
                    end
                    peak(np)=0;
                    indpeak(np)=0;
                    np=np-1;
                end
            else
                if peak(i-1)-valley<peak(i-1)/4.0
                    peak((i-1):(np-1))=peak(i:np);
                    indpeak((i-1):(np-1))=indpeak(i:np);
                    peak(np)=0;
                    indpeak(np)=0;
                    np=np-1;
                end
            end
            i=i+1;
        end
    end
    arc(1:np)=peak(1:np);
    indarc(1:np)=indpeak(1:np);
    arctheta(1:np)=theta(indarc(1:np));
    
else
    arc=0;
    indarc=0;
    arctheta=0;
end



