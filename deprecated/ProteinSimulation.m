function [PXshift] = ProteinSimulation(SmoothmRNA,Dp,Phalf,dt)
%% Calculate the amount of Protein by integrating the accumulatedmRNA (SmoothmRNA)
clear protein

%First, I need to define variables and parameters.
%gamma_P, r_p, Dp(kp)

%gamma_P : mRNA degradation rate, log(2)/half-life of mRNA (~60min)
%Phalf=60; %min, half-life of Protein, I can change this later.
gammaP=log(2)/Phalf; %Use Min. as time unit. I can change this later

%r_p:Translation rate for each AP bin, at nc14. 
rp=2; %(2 proteins/mRNA, min), I can change it later. from Petkova et al., 2014 Curr. Bio (Gregor lab)

%Dp:diffusion constant of Protein,
%Dp=0*60;    %[5:0.1:10]*60; %um^2/min. mRNA diffusion in Mouse, from BioNumbers.Yan et al., Dynamics of Translation of Single mRNA Molecules In Vivo. Cell. 2016

AP=[1:41];
dx=500/40; %delta x = 500um/40(bins) = 12.5um/bin
APbinpos = (AP-1)*dx;
kp=Dp/(dx)^2; %1/Min. unit.
% Time = mRNA.ElapsedTime;
% dt=mean(diff(Time));
%Make a zero-matrix
protein=zeros(length(SmoothmRNA),length(AP));
protein(1,:)=0;  %protein(t=0)=0 at all APbins. Initial Condition

%Additional Note.
%Think about the first and the last AP bins, as separately, 
%because they have diffusion to only one side(left or right).
%P(i,1)=P(i-1,1)+rp(1)*dt-gammaP*P(i-1,1)*dt+k*P(i-1,2)*dt-k*P(i-1,1)*dt;
%P(i,41)=...
%But, right now, we don't have to worry too much, because we don't care
%about boundary now.


%Right now, I only have data from AP bin 0.225(10th)~0.525(24th), so only use this
%datasets

%Make a Protein reaction-diffusion equation
%I will assume that the translation occurs during the whole nc14.
for i=2:length(SmoothmRNA)  %for all Timepoints
    protein(i,10)=protein(i-1,10)+rp*SmoothmRNA(i-1,10)*dt-gammaP*protein(i-1,10)*dt; %+kp*protein(i-1,11)*dt-kp*protein(i-1,10)*dt;
    protein(i,24)=protein(i-1,24)+rp*SmoothmRNA(i-1,24)*dt-gammaP*protein(i-1,24)*dt; %+kp*protein(i-1,25)*dt-kp*protein(i-1,41)*dt;
    
    for j=11:23%length(AP)-1 %for all APbins
        protein(i,j)=protein(i-1,j)+rp*SmoothmRNA(i-1,j)*dt-gammaP*protein(i-1,j)*dt+kp*dt*(protein(i-1,j-1)+protein(i-1,j+1))-2*kp*protein(i-1,j)*dt;
    end
    
end

%% Find the boundary position
APbin=[0:0.025:1];
for tpoint=3:length(SmoothmRNA)
%Spline fitting
    dxx=0.001;
    xx=0.3:dxx:0.6;
    yy=spline(APbin,protein(tpoint,:),xx);
    
    %Find the x position where it has half of the Maximum intensity
    for t=1:length(protein)
        MaxIntensity(t)=max(protein(t,:));
    end

    for i=1:length(xx)
        if yy(i)<0.5*MaxIntensity(tpoint)
            if yy(i-1)>0.5*MaxIntensity(tpoint)
                PXhalf(tpoint)=xx(i);
                HalfMaxY(tpoint)=yy(i);
            end    
        end
    end
    PXhalf(tpoint);
end

PXshift=max(PXhalf)-PXhalf(end);





end

