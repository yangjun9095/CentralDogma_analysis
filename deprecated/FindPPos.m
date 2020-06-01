function [xhalf]=FindPPos(P,t)
clear Protein
clear xx
clear yy
clear LineFit
%Analysis (mRNA)
%I need to analyze the width of the boundary
%To do that, I would need to fit a linear line, using the least square
%approximation.

%Parameters
%1. I need to pick one time point for mRNA or Protein
%For mRNA, maybe when mRNA has its maximum? and for protein as well?

%t=20;
dt=0.01;
APbin=[0:0.025:1];
%plot mRNA along the AP axis at time t.
Protein = P(round(t/dt)+1,:);





%Spline fitting
dxx=0.001;
xx=0:dxx:1;
yy=pchip(APbin,Protein,xx);
%plot(APbin,mRNA,'o',xx,yy)
%ylim([-1 400])

% for i=1:length(xx)-1
%     slope(i)=(yy(i+1)-yy(i)) / 0.001;
% end
% slope(length(xx))=slope(length(xx)-1);
% 
% Minslope=min(slope);
% pos=find(slope==Minslope);
% xpos=xx(pos);
% ypos=yy(pos);
% LineFit=Minslope*(xx-xpos)+ypos;

%Get the xpos for 0.5*Maximum
%find(yy==0.5*max(mRNA))
for i=1:length(xx)-1
    if yy(i)>0.5*max(Protein)
        if yy(i+1)<0.5*max(Protein)
            xhalf=i*dxx;
        end
    end
end
xhalf=xhalf;

% hold on
% figure(2)
% plot(xx,yy)
% plot(xhalf,yy(xhalf/dxx),'o')
% hold off
% %Get the width
% %First, get the xmax, xmin
% for i=1:length(xx)-1
%     if LineFit(i)>max(yy)
%         if LineFit(i+1)<max(yy)
%             xmax=i*dxx;
%         end
%     end
%     
% end
% 
% for i=1:length(xx)-1
%     if LineFit(i)>0
%         if LineFit(i+1)<=0|LineFit(i+1)==min(LineFit)
%             xmin=i*dxx;
%         end
%     end
%     
% end
% 
% Width = xmin-xmax;
end