function [xhalf,Width]=FindProteinPos(P,tpoint)
clear Protein
clear xx
clear yy
clear LineFit
%Analysis (Protein)
%I need to analyze the boundary position and width for Protein,'P' at specific
%time point tpoint, 'tpoint'
%To do that, I would need to fit a linear line, using the least square
%approximation.

%Parameters
%1. I need to pick one time point for Protein

APbin=[0:0.025:1];

Protein = P(tpoint,:); %Protein at specific time point


%Spline fitting
dxx=0.001;
xx=0.3:dxx:0.6;
yy=spline(APbin,Protein,xx);
%plot(APbin,mRNA,'o',xx,yy)
%ylim([-1 400])
for i=1:length(xx)-1
    slope(i)=(yy(i+1)-yy(i)) / 0.001;
end

Minslope=min(slope);
if Minslope==0
    xhalf=0;
    Width=0;
else
    pos=find(slope==Minslope);
    xpos=xx(pos);
    ypos=yy(pos);
    LineFit=Minslope*(xx-xpos)+ypos;
end

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
%Get the width
%First, get the xmax, xmin
for i=1:length(xx)-1
    if LineFit(i)>max(yy)
        if LineFit(i+1)<max(yy)
            xmax=i*dxx;
        end
    end
    
end

for i=1:length(xx)-1
    if LineFit(i)>0
        if LineFit(i+1)<=0|LineFit(i+1)==min(LineFit)
            xmin=i*dxx;
        end
    end
    
end

Width = xmin-xmax;
end