% This script is a function that calculate the mRNA profile from the given
% transcriptional rate profile (over AP, at different time points)
function [P] = CalculateProtein (M,Dp,Tp,TxlTime)
%Input parameters
%Prefix=.mat file that the mRNA profiles (over AP and time) are saved.
%M : mRNA profile (over AP bin, time)
%Dp : diffusion constant of protein (um^/sec)
%Tp : half-life of protein (min)
%TxlTime : Translational On-Time

%Assumptions:
%1) I am predicting the number of mRNA transcripts or protein molecules in
%a nucleus. Thus, I can use the MeanVectorAP (averaged over only active
%nuclei). Later, to calculate the total number of mRNA / protein over AP, I
%need to consider the fraction of active nuclei.

%2) 
%3) 

%gamma_P : protein degradation rate, log(2)/half-life of protein
gammaP=log(2)/Tp; %Use Min. as time unit.(I can change this later, for nc14)
%r_p:Tranlational rate for each AP bin, at nc14. from Petkova, 2014
rp=2; % # of protein molecules per # of mRNA molecule per min

%Dm:diffusion constant of mRNA,
Dp=Dp*60;  %um^2/min. protein diffusion %Reference : Abu-Arish, 2010 -> 7um^2/sec for Bicoid by FCS measurement

AP=1:50;
dx=500/50; %delta x = 500um/50(bins) = 10um/bin (diameter of one nucleus)
APbinpos = (AP-1)*dx+5; %Center of nuclei, defined as position of that Bin.
kp=Dp/((dx)^2); %1/Min. unit.

%Time window
dt=0.01; %min
T=[0:dt:60]; %Time matrix, spanning from 0min to 60min, with 0.01 min gap.


%Make a zero-matrix
P=zeros(length(T),length(AP));
P(1,:)=0;  %mRNA(t=0)=0 at all APbins. Initial Condition

%Additional Note.
%Think about the first and the last AP bins, as separately, 
%because they have diffusion to only one side(left or right).

%But, right now, we don't have to worry too much, because we don't care
%about boundary now.

%Make a mRNA reaction-diffusion equation
for i=2:length(T) %for all Timepoints
    for j=1:length(AP)  %for all AP bins
        
        if (T(i)<TxlTime)|(T(i)==TxlTime)
            if j==1
                P(i,1)=P(i-1,1)+rp*M(i-1,1)*dt-gammaP*P(i-1,1)*dt+kp*P(i-1,2)*dt-kp*P(i-1,1)*dt;
            elseif j==50
                P(i,50)=P(i-1,50)+rp*M(i-1,50)*dt-gammaP*P(i-1,50)*dt+kp*P(i-1,49)*dt-kp*P(i-1,50)*dt;
            else
                P(i,j)=P(i-1,j)+rp*M(i-1,j)*dt-gammaP*P(i-1,j)*dt+kp*dt*(P(i-1,j-1)+P(i-1,j+1))-2*kp*P(i-1,j)*dt;
            end
        else
            if j==1
                P(i,1)=P(i-1,1)-gammaP*P(i-1,1)*dt+kp*P(i-1,2)*dt-kp*P(i-1,1)*dt;
            elseif j==50
                P(i,50)=P(i-1,50)-gammaP*P(i-1,50)*dt+kp*P(i-1,49)*dt-kp*P(i-1,50)*dt;
            else P(i,j)=P(i-1,j)-gammaP*P(i-1,j)*dt+kp*dt*(P(i-1,j-1)+P(i-1,j+1))-2*kp*P(i-1,j)*dt;
            end
        end
    end
end

%colorful plot for timecourse : To see how the shape changes along the
%time. 
%  jStart=1;
%  jEnd=length(T)+1;
%  jTotal=jEnd-jStart;
% 
%  colormap(jet(256));
%  cmap=colormap;
%  Color=cmap(round(((jStart:jEnd)-jStart)/(jEnd-jStart)*255)+1,:);
%     
% hold on
% for l=1:500:length(T)
%     plot([0.01:0.02:0.99],M(l,:),'color',Color(l-jStart+1,:));
% %     ylim([0,1000])
% %     pause(0.01)
% end
% hold off
% 
% title('mRNA along the AP axis')
% xlabel('AP axis')
% ylabel('mRNA (molecules)')
% set(gca, 'FontSize', 30)

end