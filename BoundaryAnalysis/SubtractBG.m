%Analysis code for the CompiledParticles.mat and CompiledNuclei
%% Load the data set
%clear all

mRNA = load('CompiledParticles.mat')
Protein = load('CompiledNuclei.mat')

BG=load('BG-CompiledNuclei')
% mRNA=load('D:\Data\YangJoon\LivemRNA\Data\DynamicsResults\2016-11-23-Hb-nbGFP-MS2-NLSmCherry\CompiledParticles.mat');
% Protein=load('D:\Data\YangJoon\LivemRNA\Data\DynamicsResults\2016-11-23-Hb-nbGFP-MS2-NLSmCherry\CompiledNuclei.mat');

% mRNA2=load('D:\Data\YangJoon\LivemRNA\Data\DynamicsResults\2016-10-18-Hb-nbGFP-MS2-mCherry\CompiledParticles.mat');
% Protein2=load('D:\Data\YangJoon\LivemRNA\Data\DynamicsResults\2016-10-18-Hb-nbGFP-MS2-mCherry\CompiledNuclei.mat');

%% Subtract the nuclear fluo(w/o NB) from the nuclear fluo(w/ NB) to estimate the transcription factor concentration.

%I am assuming that the nuclear fluo without nanobody would be similar to
%the cytoplasmid fluo when the nanobody is present.

%First, I need to synchronize the two CompiledNuclei.mat
Pnc13=Protein.nc13;
Pnc14=Protein.nc14;

BGnc13=BG.nc13;
BGnc14=BG.nc14;

Tdiff=abs(Pnc13-BGnc13);
TimeDiff=BG.ElapsedTime(BG.nc13+Tdiff)-BG.ElapsedTime(BG.nc13);

%Plot from the nc13
hold on
plot(Protein.ElapsedTime(Pnc13:end)-Protein.ElapsedTime(Pnc13),Protein.MeanVectorAll(Pnc13:end),'r','LineWidth',5)
plot(BG.ElapsedTime(BGnc13:end)-TimeDiff-Protein.ElapsedTime(Pnc13),BG.MeanVectorAll(BGnc13:end),'b','LineWidth',5)

xlim([0 70])
ylim([0 150])
title('Nuclear fluorescence over time')
xlabel('Time after nc13(minute)')
ylabel('Fluorescence (AU)')
legend('Hb-Llama Tag','Control')
set(gca,'Fontsize',30)
%% Plot only during the nc14
hold on
plot(Protein.ElapsedTime(Pnc14:end)-Protein.ElapsedTime(Pnc14),Protein.MeanVectorAll(Pnc14:end),'r','LineWidth',5)
plot(BG.ElapsedTime(BGnc14:end)-TimeDiff-Protein.ElapsedTime(Pnc14),BG.MeanVectorAll(BGnc14:end),'b','LineWidth',5)

%xlim([0 70])
ylim([0 150])
title('Nuclear fluorescence over time')
xlabel('Time after nc14(minute)')
ylabel('Fluorescence (AU)')
legend('Hb-Llama Tag','Control')
set(gca,'Fontsize',30)

%% Subtracting the nuc.fluo of control (Time synch + same AP bin)

%The length of measurement can be different.
%I will choose the shorter one as reference.
%Compare the length of Timeseries in nc14.
if length(Protein.ElapsedTime(Pnc14:end)) > length(BG.ElapsedTime(BGnc14:end))
    for i=1:length(BG.ElapsedTime(BGnc14:end))
        NewProtein(i,:)=Protein.MeanVectorAP(Pnc14+i+Tdiff-1,:)-BG.MeanVectorAP(BGnc14+i-1,:);
    end
else 
    for i=1:length(Protein.ElapsedTime(BGnc14:end))
        NewProtein(i,:)=BG.MeanVectorAP(BGnc14+i-1,:)-Protein.MeanVectorAP(Pnc14+i+Tdiff-1,:);
    end
end

hold on
%plot(Protein.ElapsedTime(Pnc14:end)-Protein.ElapsedTime(Pnc14),Protein.MeanVectorAll(Pnc14:end),'r','LineWidth',5)
plot(BG.ElapsedTime(BGnc14:end)-TimeDiff-BG.ElapsedTime(BGnc14),NewProtein,'LineWidth',5)

%xlim([0 70])
%ylim([0 150])
title('Nuclear fluorescence over time')
xlabel('Time after nc14(minute)')
ylabel('Fluorescence (AU)')
%legend('Hb-Llama Tag','Control')
set(gca,'Fontsize',30)

%% NewProtein -> Analyze how (1)the boundary position and (2) the boundary sharpness change over time.