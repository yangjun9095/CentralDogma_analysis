function [Xhalf]=FindPos(M,tpoint)
clear mRNA
clear xx
clear yy
%clear LineFit
%Analysis (mRNA for Puresimulaton data)
%I need to analyze the width of the boundary
%To do that, I would need to fit a linear line, using the least square
%approximation.

%Parameters
%1. I need to pick one time point for mRNA or Protein
%For mRNA, maybe when mRNA has its maximum? and for protein as well?

%t=20; %min
%dt=0.01;
APbin=[0:0.025:1];
%plot mRNA along the AP axis at time t.
mRNA = M(tpoint,:);

%Spline fitting
dxx=0.001;
xx=0.2:dxx:0.6;
yy=spline(APbin,mRNA,xx);
plot(APbin,mRNA,'o',xx,yy)
xlim([0.2 0.6])

% for i=1:length(xx)-1
%     slope(i)=(yy(i+1)-yy(i)) / 0.001;
% end

% Minslope=min(slope);
% pos=find(slope==Minslope);
% xpos=xx(pos);
% ypos=yy(pos);
% LineFit=Minslope*(xx-xpos)+ypos;

%Get the xpos for 0.5*Maximum
%find(yy==0.5*max(mRNA))
MaxIntensity(tpoint)=max(mRNA);

for i=1:length(xx)
    if yy(i)<0.5*MaxIntensity(tpoint)
        if yy(i-1)>0.5*MaxIntensity(tpoint)
            Xhalf=xx(i);
        end    
    end
end
%Get the width
%First, get the xmax, xmin
% for i=1:length(xx)-1
%     if LineFit(i)>max(yy)
%         if LineFit(i+1)<max(yy)
%             xmax=i*dxx;
%         end
%     end
%     
% end

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