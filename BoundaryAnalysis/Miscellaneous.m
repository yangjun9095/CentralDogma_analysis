%Miscellaneous
%% Plotting spot fluorescence for time (for different AP position as different color)
%Use MeanVectorAP at APbinID = ap, and plot them for different AP bins
jStart=1; %Protein.nc14;
jEnd=41;
jTotal=jEnd-jStart;

colormap(jet(64));
cmap=colormap;
Color=cmap(round(((jStart:jEnd)-jStart)/(jEnd-jStart)*63)+1,:);

hold on
for ap=1:41
    plot(mRNA.ElapsedTime,mRNA.MeanVectorAP(:,ap),'color',Color(ap-jStart+1,:))
    drawnow
    pause(1)
end
%plot(mRNA.ElapsedTime(mRNA.nc13:end)-(mRNA.nc13-mRNA.nc13),mRNA.MeanVectorAll(mRNA.nc13:end),'b')
hold off

title('Transcription Rate over time')
xlabel('Time(min)')
ylabel('Transcription Rate')