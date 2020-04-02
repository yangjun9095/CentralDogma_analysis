% This script is a function that calculate the mRNA profile from the given
% transcriptional rate profile (over AP, at different time points)
function [M] = CalculatemRNA (Rate,Dm,Tm,TxnTime)
%Input parameters
%Prefix=.mat file that the mRNA profiles (over AP and time) are saved.
%Rate : Rate of transcription initiation (over AP bin, time)
%Dm : diffusion constant of mRNA (um^/sec)
%Tm : half-life of mRNA (min)
%TxnTime : Transcriptional On-Time

%Assumptions:
%1) I am predicting the number of mRNA transcripts or protein molecules in
%a nucleus. Thus, I can use the MeanVectorAP (averaged over only active
%nuclei). Later, to calculate the total number of mRNA / protein over AP, I
%need to consider the fraction of active nuclei.

%2) I will also assume that the transcription will persist for only 20min
%after the onset of nc 14.

%3) 

%gamma_M : mRNA degradation rate, log(2)/half-life of mRNA (~60min)
gammaM=log(2)/Tm; %Use Min. as time unit.(I can change this later, for nc14)
%r_m:Transcription rate for each AP bin, at nc14. from Hernan's 2013 Curr.Bio paper
rm=Rate; %Normalized by # of RNAP

%Dm:diffusion constant of mRNA,
Dm=Dm*60;  %um^2/min. mRNA diffusion in Mouse, from BioNumbers.Yan et al., Dynamics of Translation of Single mRNA Molecules In Vivo. Cell. 2016

AP=1:50;
dx=500/50; %delta x = 500um/50(bins) = 10um/bin (diameter of one nucleus)
APbinpos = (AP-1)*dx+5; %Center of nuclei, defined as position of that Bin.
km=Dm/((dx)^2); %1/Min. unit.

%Time window
dt=0.01; %min
T=[0:dt:60]; %Time matrix, spanning from 0min to 60min, with 0.01 min gap.

%We need Transcription Time Window, along the AP axis, as in Fig.4B (1x41)
%Twindow=[0,0,0,0,0,0,0,0,18,17,17.5,17,19,17,15.5,16,15.5,12.5,11,11.5,10,8.5,11,8.5,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6];

% %Fraction of Active Nuclei Fig.4.F
% FractionActiveNuclei=MaxOnRatio(:,3);
% FractionActiveNuclei(isnan(FractionActiveNuclei))=0;

%Effective rate of transcription along the AP axis.
%rmeff=rm.*FractionActiveNuclei;

%Make a zero-matrix
M=zeros(length(T),length(AP));
M(1,:)=0;  %mRNA(t=0)=0 at all APbins. Initial Condition

%Additional Note.
%Think about the first and the last AP bins, as separately, 
%because they have diffusion to only one side(left or right).

%But, right now, we don't have to worry too much, because we don't care
%about boundary now.

%Make a mRNA reaction-diffusion equation
for i=2:length(T) %for all Timepoints
    for j=1:length(AP)  %for all AP bins
        
        if T(i)<TxnTime
            if j==1
                M(i,1)=M(i-1,1)+rm(j)*dt-gammaM*M(i-1,1)*dt+km*M(i-1,2)*dt-km*M(i-1,1)*dt;
            elseif j==50
                M(i,50)=M(i-1,50)+rm(j)*dt-gammaM*M(i-1,50)*dt+km*M(i-1,49)*dt-km*M(i-1,50)*dt;
            else
            M(i,j)=M(i-1,j)+rm(j)*dt-gammaM*M(i-1,j)*dt+km*dt*(M(i-1,j-1)+M(i-1,j+1))-2*km*M(i-1,j)*dt;
            end
        else
            if j==1
                M(i,1)=M(i-1,1)-gammaM*M(i-1,1)*dt+km*M(i-1,2)*dt-km*M(i-1,1)*dt;
            elseif j==50
                M(i,50)=M(i-1,50)-gammaM*M(i-1,50)*dt+km*M(i-1,49)*dt-km*M(i-1,50)*dt;
            else M(i,j)=M(i-1,j)-gammaM*M(i-1,j)*dt+km*dt*(M(i-1,j-1)+M(i-1,j+1))-2*km*M(i-1,j)*dt;
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