% This script is a function that calculate the mRNA profile from a given
% transcriptional rate profile (over AP, at different time points)
function [AccumulatedmRNA,Time] = estimatemRNA (Rate, Dm, Tm, Time, Tmax)
%Input parameters
% Prefix=.mat file that the mRNA profiles (over AP and time) are saved.
% Rate : Rate of transcription initiation (over AP bin, time)
% Dm : diffusion constant of mRNA (um^/sec)
% Tm : half-life of mRNA (min)
% Time : ElapsedTime corresponding to the Rate for interpolation (we need to
% plug in as the unit of [minutes].
% Tmax : maximum time length

%Assumptions:
%1) I am predicting the number of mRNA transcripts or protein molecules @ each AP bin. 
% Thus, I can use the MeanVectorAP (averaged over only active nuclei) * 
% NParticlesAP (number of particles detected). 
% This means that I need to plug in "Rate" as "MeanVectorAP.*NParticlesAP"

% Define parameters
%gamma_M : mRNA degradation rate, log(2)/half-life of mRNA (~60min)
gammaM=log(2)/Tm; %Use Min. as time unit.(I can change this later, for nc14)
%r_m:Transcription rate for each AP bin, at nc14. from Hernan's 2013 Curr.Bio paper
%rm=Rate; %Normalized by # of RNAP

%Dm:diffusion constant of mRNA,
Dm=Dm*60;  %um^2/min. 
%mRNA diffusion in Mouse, from BioNumbers.Yan et al., Dynamics of Translation of Single mRNA Molecules In Vivo. Cell. 2016

AP=1:40;
dx=500/40; %delta x = 500um/50(bins) = 10um/bin (diameter of one nucleus)
APbinpos = (AP-1)*dx + 1/2*dx; %Center of nuclei, defined as position of that Bin.
km=Dm/((dx)^2); %1/Min. unit.

% Time window
% The condition to do a numerical calculation is that
% dt << 1/k, thus dt should be set by the upper limit of k (which is
% basically the upper limit of Dm). In here, I'll just assume the maximum
% Dm as 100um^2/sec (as this is crazy fast). Then, the 1/k = 0.026, so dt
% should be smaller than this.
dt=0.01; %min
T=[0:dt:Tmax]; %Time matrix, spanning from 0min to Tmax(i.e.60min), with dt (min) interval.

% Interpolate the "rm" from the given input.
rm = interp1(Time, Rate, T); % rm is now interpolated with dt.

%Make a zero-matrix
M=zeros(length(T),length(AP));
M(1,:)=0;  %mRNA(t=0)=0 at all APbins. Initial Condition


%Additional Note.
%Think about the first and the last AP bins, as separately, 
%because they have diffusion to only one side(left or right).
%But, right now, we don't have to worry too much, because we don't care
%about boundary now.

%Write down a chemical master equation of a reaction-diffusion system
for i=2:length(T) %for all Timepoints
    for j=1:length(AP)  %for all AP bins       
        if j==1
            M(i,1)=M(i-1,1)+rm(j)*dt-gammaM*M(i-1,1)*dt+km*M(i-1,2)*dt-km*M(i-1,1)*dt;
        elseif j==max(AP)
            M(i,max(AP))=M(i-1,max(AP))+rm(j)*dt-gammaM*M(i-1,max(AP))*dt+km*M(i-1,max(AP)-1)*dt-km*M(i-1,max(AP))*dt;
        else
        M(i,j)=M(i-1,j)+rm(j)*dt-gammaM*M(i-1,j)*dt+km*dt*(M(i-1,j-1)+M(i-1,j+1))-2*km*M(i-1,j)*dt;
        end
    end
end

%% Deprecated
%We need Transcription Time Window, along the AP axis, as in Fig.4B (1x41)
%in Garcia, 2013?
%Twindow=[0,0,0,0,0,0,0,0,18,17,17.5,17,19,17,15.5,16,15.5,12.5,11,11.5,10,8.5,11,8.5,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6];

% %Fraction of Active Nuclei Fig.4.F
% FractionActiveNuclei=MaxOnRatio(:,3);
% FractionActiveNuclei(isnan(FractionActiveNuclei))=0;

%Effective rate of transcription along the AP axis.
%rmeff=rm.*FractionActiveNuclei;



%% colorful plot for the time-course : To see how the shape changes along the time. 
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