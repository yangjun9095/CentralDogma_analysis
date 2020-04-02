


for j=1:length(T)-1
    [xhalf,Width]=FindPos(M,j);
    Result(j,1)=xhalf;
    Result(j,2)=Width;
end
Result