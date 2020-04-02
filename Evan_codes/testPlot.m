function testPlot(time)

tempInterMeasProteinFile = load('E:\EvanM\RandomVariables\tempInterMeasProtein','tempInterMeasProtein');
tempInterMeasProtein = tempInterMeasProteinFile.tempInterMeasProtein;

y = tempInterMeasProtein(time,:);

figure
hold on
plot(y)
hold off