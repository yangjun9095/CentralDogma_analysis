function [xhalf,Width,Slope]=FindPosition(Profile,tpoint)
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

Profile = Profile(tpoint,:); %Protein(or mRNA) at specific time point

% Smoothen the profile over 3 AP windows, neighboring AP bins (7.5%)
NewProfile = nan(1,41);
for i=9:33 % 20%~80% AP axis
    NewProfile(i) = nanmean(Profile(i-1:i+1));
end

%pchip interpolation
dxx=0.001;
xx=0.2:dxx:0.6;
yy=pchip(APbin(9:25),NewProfile(9:25),xx);

% Plot to check if the pchip fitting is okay
%plot(APbin,Profile,'o',xx,yy)
%ylim([-1 400])
for i=21:length(xx)-21
    slope(i)=(yy(i+20)-yy(i-20)) / (40*dxx);
end

Minslope=min(slope);
% 
% if Minslope==0
%     xhalf=0;
%     Width=0;
%     Slope = 0;
% else
%     pos=find(slope==Minslope);
%     xpos=xx(pos);
%     ypos=yy(pos);
%     LineFit=Minslope*(xx-xpos)+ypos;
% end

%Get the xpos for 0.5*Maximum
%find(yy==0.5*max(mRNA))
if Minslope ~= 0 
    for i=1:length(xx)-1
        if yy(i)>0.5*max(NewProfile(9:25)) && yy(i+1)<0.5*max(NewProfile(9:25))
            xhalf=xx(i);
            Slope = slope(i);
            %     pos=find(slope==Minslope);
            xpos=xx(i);
            ypos=yy(i);
%             xhalf=xhalf;
%             Slope =Slope;
%         else
%             error('there is no boundary? check the profile')
        end
    end
    LineFit=Slope*(xx-xpos)+ypos;
    %Get the width
    %First, get the xmax, xmin
    for i=1:length(xx)-1
        if LineFit(i)>=max(yy) && LineFit(i+1)<=max(yy)
            xmax=xx(i);
        end
    end

    for i=1:length(xx)-1
        if LineFit(i)>=0 && (LineFit(i+1)<=0 || LineFit(i+1)==min(LineFit))
            xmin=xx(i);
        end
    
    end
    Width = xmin-xmax;
    
else
    xhalf=0.5;
    Width=1;
    Slope = 0;
    
end


end