function [M]=SimulatemRNA(Dm)

load('DataForPlots.mat');
clear M
%Rate = RateWeight;
for i=1:41
    if i<18
        Rate(i)=7000;
    elseif i>18
        Rate(i)=0;
    else Rate(i)=3500;
    end
end

%Let's start with considering nc14 only.
%First, I need to define variables and parameters.

%gamma_M : mRNA degradation rate, log(2)/half-life of mRNA (~60min)
gammaM=log(2)/60; %Use Min. as time unit.
%r_m:Transcription rate for each AP bin, at nc14. from Hernan's 2013 Curr.Bio paper
rm=Rate; %Normalized by # of RNAP
%Remove Nans
for i=1:41
    rm(isnan(rm))=0;
end
%Dm:diffusion constant of mRNA,
%Dm=0.1*60; %[0.1:0.1:0.4]*60; %um^2/min. mRNA diffusion in Mouse, from BioNumbers.Yan et al., Dynamics of Translation of Single mRNA Molecules In Vivo. Cell. 2016

AP=[1:41];
APbin=[0:0.025:1];
dx=500/40; %delta x = 500um/40(bins) = 12.5um/bin
APbinpos = (AP-1)*dx;
k=Dm/(dx)^2; %1/Min. unit.

%Time window
dt=0.01;
T=[0:dt:60]; %Time matrix, spanning from 0min to 60min, with 0.01 min gap.
%We need Transcription Time Window, along the AP axis, as in Fig.4B (1x41)
Twindow=[0,0,0,0,0,0,0,0,18,17,17.5,17,19,17,15.5,16,15.5,12.5,11,11.5,10,8.5,11,8.5,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6];

%Fraction of Active Nuclei Fig.4.F
FractionActiveNuclei=MaxOnRatio(:,3);
FractionActiveNuclei(isnan(FractionActiveNuclei))=0;

%Effective rate of transcription along the AP axis. (1x41)
rmeff=transpose(rm).*FractionActiveNuclei;

%Make a zero-matrix
M=zeros(length(T),length(AP));
M(1,:)=0;  %mRNA(t=0)=0 at all APbins. Initial Condition

%Make a mRNA reaction-diffusion equation
for j=2:length(AP)-1 %for all APbins
    for i=2:length(T)  %for all Timepoints
        if i<Twindow(j)
            M(i,j)=M(i-1,j)+rm(j)*dt-gammaM*M(i-1,j)*dt+k*dt*(M(i-1,j-1)+M(i-1,j+1))-2*k*M(i-1,j)*dt;
        else
            M(i,j)=M(i-1,j)-gammaM*M(i-1,j)*dt+k*dt*(M(i-1,j-1)+M(i-1,j+1))-2*k*M(i-1,j)*dt;
        end
    end
end



end