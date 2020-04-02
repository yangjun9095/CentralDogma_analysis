function [xpos] = FindProteinPosition(Protein,t)
%I need to analyze the width of the boundary
%To do that, I would need to fit a linear line, using the least square
%approximation.

%Parameters
%1. I need to pick one time point for mRNA or Protein
%For mRNA, maybe when mRNA has its maximum? and for protein as well?

%t=100 %Time Index

%plot Protein along the AP axis at time t.
%Protein = SmoothProtein(t,:);
%MAX = max(mRNA);
%mRNA = mRNA/MAX; %for normalization



% hold on
%plot ([0:0.025:1],mRNA,'o') % Normalized mRNA along the AP axis
%plot ([0:0.025:1-0.025],slope,'o') %slope
% hold off
APbin=[0:0.025:1];
%hold on
%Spline fitting
xx=0.2:0.001:0.6;
yy=spline(APbin,Protein,xx);
% plot(APbin,Protein,'o',xx,yy)
% ylim([-1 2000])
for i=1:length(xx)-1
    slope(i)=(yy(i+1)-yy(i)) / 0.001;
end

%Use the half of the maximum point's x value's slope as a point where I set
%position!

pos = sum(yy>=0.5*max(Protein));
xpos=xx(pos);
ypos=yy(pos);
%plot(xx,LineFit)

% title('Protein along AP axis','fontsize',25)
% xlabel('AP axis','fontsize',20)
% ylabel('Protein concentration(AU)','fontsize',20)
% legend('Protein')

% barh(max(yy),1);
% barh(0,1);

hold off
end