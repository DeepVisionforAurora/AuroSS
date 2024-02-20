% 绘制日侧0600-1800MLT时间范围内的AACGM坐标系
% 纬度范围：[MinLat 90] 纬度间隔：5
% MinLat：最小纬度值
% Ver:2006.06.06.002
%----------------------------------------------

function AACGM_figure=half_Axis4AACGM()
MinLat=65;
L=MinLat:5:90;
r=90-L;
[r_m,r_n]=size(r);
Theta=0:pi/360:pi;

% AACGM_figure=figure;
% set(AACGM_figure,'Position',[360 33 650 650]);
% AACGM_axis=get(AACGM_figure,'Children');
% set(AACGM_axis,'Position',[0.11 0.11 0.815 0.815]);
% set(AACGM_axis,'XLim',[-27 27],'YLim',[-27 27]);
% set(AACGM_axis,'XTick',[],'XTickLabel',[]);
% set(AACGM_axis,'YTick',[],'YTickLabel',[]);
%----------------------------------------------
% 绘制纬度圈
for ii=1:r_n
    [x,y]=pol2cart(Theta,cos(0*Theta)*r(1,ii));
    plot(x,y,'-k');
    hold on;
    clear x y;
end
%----------------------------------------------
MLT=0:pi/12:pi;
MLTStr=['06';'07';'08';'09';'10';'11';'12';'13';'14';'15';'16';'17';'18'];
%           6     8       10      12    14    16      18
MLTStr_Cof=[0 1/6 2/6 3/6 4/6 5/6 1 4/3 5/3 2 7/3 8/3 3 ];

[MLT_m,MLT_n]=size(MLT);
MLat=90-[MinLat:0.01:90];
[MLat_m,MLat_n]=size(MLat);
E=ones(MLat_m,MLat_n);
% 绘制MLT分界线
for jj=1:MLT_n
    [x,y]=pol2cart(MLT(1,jj)*E,MLat);
%     [tx,ty]=pol2cart(MLT(1,jj)*E,90-(60-MLTStr_Cof(jj)));
    plot(x,y,'-k');
%     [x_m,x_n]=size(x);
%     text(tx,ty,MLTStr(jj,:),'FontSize',12); % MLT时间标识
    hold on
    clear x y;
end
%----------------------------------------------
text(-21,-1,'70\circ','FontName','Times New Roman','FontSize',12,'FontWeight','bold');
text(-11,-1,'80\circ','FontName','Times New Roman','FontSize',12,'FontWeight','bold');
text(-1,-1,'90\circ','FontName','Times New Roman','FontSize',12,'FontWeight','bold');
text(9,-1,'80\circ','FontName','Times New Roman','FontSize',12,'FontWeight','bold');
text(19,-1,'70\circ','FontName','Times New Roman','FontSize',12,'FontWeight','bold');

text(-27,1,'18','FontName','Times New Roman','FontSize',10,'FontWeight','bold');
text(25.5,1,'06','FontName','Times New Roman','FontSize',10,'FontWeight','bold');
text(-0.7,26,'12','FontName','Times New Roman','FontSize',10,'FontWeight','bold');

% text(5,37,'AACGM Coordinate','FontSize',12,'FontWeight','Bold');
axis equal tight;
AACGM_axis=get(gcf,'Children');
% % set(AACGM_axis,'Position',[0.11 0.11 0.815 0.815]);
set(gca,'XLim',[-28 28],'YLim',[-3 28]);
set(gca,'XTick',[],'XTickLabel',[]);
set(gca,'YTick',[],'YTickLabel',[]);