function [Xhalf,Width] = GetBoundary_Dostatni (Profile)
%Description : 
% This is a function to calculate the boundary position and width over
% time, as in Crauk and Dostatni, 2005 paper. Detailed description could be
% found in the paper.

%Input paramters
% Profile : Profile or protein level over AP, rows-time, column-AP bins
%

%% To capture the boundary position and width at one specific time point

%I need to interpolate the x axis (AP), and I will use pchip for
%interpolation. Alternatives : interp1, 

%Prifle=cell2mat(Profile(5)); Profile should be a matrix, if it's an element of a cell,
%it should be converted by cell2mat.

%hold on
APbin=0.01:0.02:0.99;
for T=2:length(Profile(1,:)) %6001
    Profile(T,:)=Profile(T,:);
    %plot([0.01:0.02:0.99],Profile(T,:))
    
    clear Left
    clear Right
    
    clear yyy
    %At one time point, T, get the boundary position and width
    
    %First, I need to do interpolation(for space) using pchip
    xxx=0:0.001:1; % Interpolate with 0.001 (0.1% of APbin)
    
    yyy(T,:)=pchip(APbin,Profile(T,:),xxx);
    
    for i=1:length(xxx)
        if (yyy(T,i)<0.5*max(Profile(T,:)))|(yyy(T,i)==0.5*max(Profile(T,:)))
            if yyy(T,i-1)>0.5*max(Profile(T,:))
                Xhalf(T)=xxx(i);
                Yhalf(T)=yyy(T,i);
                Slope(T)=(yyy(T,i+2)-yyy(T,i-2))/(xxx(i+2)-xxx(i-2));
            end
        else
            Xhalf(T)=0;
            Yhalf(T)=0;
            Slope(T)=0;
        end
        
        %Tangential line at the boundary position
        yy(T,:)=Slope(T)*(xxx-Xhalf(T))+Yhalf(T);
    end
    
    for j=1:length(xxx)
        if (yy(T,:)<max(Profile(T,:)))+(yy(T,:)==max(Profile(T,:))) %In case the tangent line is always smaller than the Profile profile
            Left=0;
        elseif yy(T,j)>max(Profile(T,:))
            if yy(T,j+1)<max(Profile(T,:))
                Left=j;
            end
        end
        
        if (yy(T,:)>min(Profile(T,:)))+(yy(T,:)==min(Profile(T,:)))
            Right=0;
        elseif yy(T,j)<min(Profile(T,:))
            if yy(T,j-1)>min(Profile(T,:))
                Right=j;
            end
        end
        
    end
    
    Width(T)=(Right-Left)*0.001;
    Sharp(T)=abs(Slope(T));
    
end
% %hold off
% %% Plot to check Xhalf and Width
% hold on
% plot(0:0.01:60,Xhalf)
% plot(0:0.01:60,Width)
% %plot(0:0.01:60,Width2)
% %% Plot to check curve fitting
% hold on
% plot(xxx,yy(5757,:))
% plot(xxx,yyy(5757,:))



end