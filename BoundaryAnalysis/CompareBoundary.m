%Comparison of Boundary position (Xhalf, PXhalf(for different parameters), and NBXhalf) ,
%and Width (Width, PWidth(for diff. parameters), and NBWidth)

%First, I will plot Boundary Position for AccumulatedmRNA, Protein (with
%different parameters, 

PDiffusion = [0,1,5,10];
Phalf=[3,6,60,600];

%% Simulation for different Dp
[PXhalf1,PWidth1] = PredictProtein(PDiffusion(1),Phalf(1),NewTime,InterpSmoothmRNA);
[PXhalf2,PWidth2] = PredictProtein(PDiffusion(2),Phalf(1),NewTime,InterpSmoothmRNA);
[PXhalf3,PWidth3] = PredictProtein(PDiffusion(3),Phalf(1),NewTime,InterpSmoothmRNA);
[PXhalf4,PWidth4] = PredictProtein(PDiffusion(4),Phalf(1),NewTime,InterpSmoothmRNA);

%% Simulation for different P half-life
[PXhalf5,PWidth5] = PredictProtein(PDiffusion(3),Phalf(1),NewTime,InterpSmoothmRNA);
[PXhalf6,PWidth6] = PredictProtein(PDiffusion(3),Phalf(2),NewTime,InterpSmoothmRNA);
[PXhalf7,PWidth7] = PredictProtein(PDiffusion(3),Phalf(3),NewTime,InterpSmoothmRNA);
[PXhalf8,PWidth8] = PredictProtein(PDiffusion(3),Phalf(4),NewTime,InterpSmoothmRNA);
%% Plot for Boundary position - different Dp
hold on
plot(SmoothTime,Xhalf,'r')
plot(NewTime,PXhalf1)
plot(NewTime,PXhalf2)
plot(NewTime,PXhalf3)
plot(NewTime,PXhalf4)
plot(NBTime(40:end),NBXhalf(39:end),'o')
hold off
title('Comparison of Boundary position','fontsize',30)
xlabel('Time (min)')
ylabel('AP')
ylim([0.3 0.5])
legend('AccumulatedmRNA','Dp=0\mum ^{2}/sec','Dp=1\mum ^{2}/sec','Dp=5\mum ^{2}/sec','Dp=10\mum ^{2}/sec','Nanobody Data')
set(gca,'fontsize',25)

%% Plot for Boundary position - different P half
hold on
plot(SmoothTime,Xhalf,'r')
plot(NewTime,PXhalf5)
plot(NewTime,PXhalf6)
plot(NewTime,PXhalf7)
plot(NewTime,PXhalf8)
plot(NBTime(40:end),NBXhalf(39:end),'o')
hold off
title('Comparison of Boundary position','fontsize',30)
xlabel('Time (min)')
ylabel('AP')
ylim([0.3 0.5])
legend('AccumulatedmRNA','Half-life=3min','Half-life=6min','Half-life=60min','Half-life=600min','Nanobody Data')
set(gca,'fontsize',25)


%% Plot for Boundary Width - Dp
hold on
plot(SmoothTime,Width,'r')
plot(NewTime,PWidth1)
plot(NewTime,PWidth2)
plot(NewTime,PWidth3)
plot(NewTime,PWidth4)
plot(NBTime(40:end),NBWidth(37:end),'o')
hold off
title('Comparison of Boundary Width','fontsize',30)
xlabel('Time (min)')
ylabel('AP')
%ylim([0.3 0.5])
legend('AccumulatedmRNA','Dp=0\mum ^{2}/sec','Dp=1\mum ^{2}/sec','Dp=5\mum ^{2}/sec','Dp=10\mum ^{2}/sec','Nanobody Data','Location','northwest')
set(gca,'fontsize',25)

%% Plot for Boundary Width - Half-life
hold on
plot(SmoothTime,Width,'r')
plot(NewTime,PWidth5)
plot(NewTime,PWidth6)
plot(NewTime,PWidth7)
plot(NewTime,PWidth8)
plot(NBTime(40:end),NBWidth(37:end),'o')
hold off
title('Comparison of Boundary Width','fontsize',30)
xlabel('Time (min)')
ylabel('AP')
%ylim([0.3 0.5])
legend('AccumulatedmRNA','Half-life=3min','Half-life=6min','Half-life=60min','Half-life=600min','Nanobody Data','Location','northwest')
set(gca,'fontsize',25)