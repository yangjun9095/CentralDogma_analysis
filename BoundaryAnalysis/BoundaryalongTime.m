%% Find the mRNA boundary position & Width along the timecourse

hold on
for i=2:length(T)-1
    xhalf(i)=FindPPos(M,T(i));
    
    %plot(T(i),Width,'o')
end
hold off
plot(T(2:end),xhalf,'o','color','r')

title('Boundary Width vs Time')
xlabel('Time(min)')
ylabel('Boundary Width (AP)')
ylim([0,1])
    
% title('Boundary Position vs Time')
% xlabel('Time(min)')
% ylabel('Boundary Position (AP)')
% ylim([0,1])

%% Find the Protein boundary position & Width along the timecourse
xhalf=zeros(1,length(T)-1);

for i=10:length(T)-1
    xhalf(i)=FindPPos(P,T(i));    
end

plot(T(10:end-1),xhalf(10:end)*50,'o')
    %plot(T(i),Width,'o')
% title('Boundary Width vs Time')
% xlabel('Time(min)')
% ylabel('Boundary Width (AP)')
%ylim([0,1])
    
title('Boundary Position vs Time')
xlabel('Time(min)')
ylabel('Boundary Position (# of Nuclei)')
%ylim([0,1])
legend('Simulation')
set(gca,'fontsize',20)