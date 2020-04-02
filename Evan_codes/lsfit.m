function [difference] = lsfit(measuredData,predictedData,t)

digits(8)

difference = 0;
if size(measuredData,1) == size(predictedData,1) && size(measuredData,2) == size(predictedData,2)
    for j = 1:size(measuredData,2)
            difference = difference + (measuredData(t,j) - predictedData(t,j))^2;
    end
else
    error('measuredData and predictedData are not the same size')
end
        