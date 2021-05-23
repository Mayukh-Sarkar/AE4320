function [pks,pks_locs] = PeakLabels(x,y,hs,vs)

[pks,pks_locs] = findpeaks(y,x); % find peaks and their locations
figure('WindowState','maximized') % start a maximized figure window
plot(x,y,'LineWidth',2,'Color','b'); % plot input arrays
hold on;
plot(pks_locs,pks,'ok','LineWidth',1.5,'MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',10) % plot peaks in their locations
hold on
ylim([min(y)-0.5 max(y)+0.5]) % leave some y margin
text(pks_locs+hs,pks+vs,num2str(pks_locs','%.1f'),'FontSize',16); % place dip labels
box on
% lines below control some plot aspects to give them a nice look
set(gca, 'Layer', 'top')
set(gca,'fontsize',28) 
set(gca,'linewidth',3)
set(gca,'FontName','Times New Roman')
set(gcf, 'Position', [60,40,1150,600])
set(gcf,'color','w'); %%%
pbaspect([2 1 1])

end