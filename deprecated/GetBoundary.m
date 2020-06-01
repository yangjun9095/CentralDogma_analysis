%This is a function to calculate the boundary position and width over time.
function [Xhalf,Width,yyy,Left,Right,Yhalf,yy] = GetBoundary (M)

%Input paramters
%M : mRNA profile, rows-time, column-AP bins -> Although this is for mRNA,
%this can also be applied for protein
%

%Note. For now, I will get the boundary for 5min, 10min,... as 5min
%interval

%% To capture the boundary position and width at one specific time point

%I need to interpolate the x axis (AP), and I will use pchip for
%interpolation. Alternatives : interp1, 

%M=cell2mat(mRNA(5)); M should be a matrix, if it's an element of a cell,
%it should be converted by cell2mat.
clear Left
clear Right
Left=[];
Right=[];
clear yyy
yyy=[];

%hold on
APbin=0.01:0.02:0.99;
for T=501:500:6001
    MRNA(T,:)=M(T,:);
    %plot([0.01:0.02:0.99],MRNA(T,:))

    %At one time point, T, get the boundary position and width
    
    %First, I need to do interpolation(for space) using pchip
    xxx=0:0.001:1;
    
    yyy(T,:)=pchip(APbin,MRNA(T,:),xxx);
    
    for i=1:length(xxx)
        if (yyy(T,i)<0.5*max(MRNA(T,:)))|(yyy(T,i)==0.5*max(MRNA(T,:)))
            if yyy(T,i-1)>0.5*max(MRNA(T,:))
                Xhalf(T)=xxx(i);
                Yhalf(T)=yyy(T,i);
                Slope(T)=(yyy(T,i+1)-yyy(T,i-1))/(xxx(i+1)-xxx(i-1));
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
        if (yy(T,:)<max(MRNA(T,:)))+(yy(T,:)==max(MRNA(T,:))) %In case the tangent line is always smaller than the mRNA profile
            Left(T)=0;
        elseif yy(T,j)>max(MRNA(T,:))
            if yy(T,j+1)<max(MRNA(T,:))
                Left(T)=j;
            end
        end
        
        if (yy(T,:)>min(MRNA(T,:)))+(yy(T,:)==min(MRNA(T,:)))
            Right(T)=0;
        elseif yy(T,j)<min(MRNA(T,:))
            if yy(T,j-1)>min(MRNA(T,:))
                Right(T)=j;
            end
        end
        
    end
    
    Width(T)=(Right(T)-Left(T))*0.001;
    Sharp(T)=abs(Slope(T));
    
%     hold on
%     plot(0.01:0.02:0.99,M(T,:),'o')
%     plot(xxx,yyy(T,:))
    
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