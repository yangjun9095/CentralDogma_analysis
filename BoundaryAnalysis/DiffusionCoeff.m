Dm=[0.1:0.1:1];
hold on
for i=1:length(Dm)
    mRNA = SimulatemRNA(Dm(i));
    mRNA=mRNA(600,:);
    plot([0:0.025:1],mRNA)
end
hold off
    