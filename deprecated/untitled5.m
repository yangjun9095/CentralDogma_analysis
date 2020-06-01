clear ResultP

Dp=[1:1:10]*60; %um^2/sec

for i=1:length(Dp)
PROTEIN(i)={SimulateProtein(Dp(i))};
end

for j=1:length(Dp)
    clear Protein
    [xhalf,Width]=FindPPos(cell2mat(PROTEIN(j)),20);
    ResultP(j,1)=xhalf;
    ResultP(j,2)=Width;
end
ResultP
