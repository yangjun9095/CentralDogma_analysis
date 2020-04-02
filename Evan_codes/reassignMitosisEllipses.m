%%Reassigns the positions of the nuclei during the first half of mitosis to
%%the last clear positional assignment. Does the same for the second half
%%by reassigning them to the first clear image in the next nuclear cycle.

%%Convention is that if the number of mitotic frames is odd, the second
%%half i.e. the half being assigned to the next nuclear cycle, is longer.

%variable 'first' is first bad frame and variable 'last' is last bad frame.
%'filePath' variable is Prefix



function reassignMitosisEllipses(filePath,first,last)

%%determine file path
filePath1 = 'E:\EvanM\LivemRNA\Data\DynamicsResults\';
filePath2 = '\Ellipses.mat';
s = strcat(filePath1,filePath,filePath2);

interval = last-first;
roundedHalfWay = round((interval)/2);

%Load Ellipses
load(s)

%Replace Ellipses
for i=first:first+roundedHalfWay-1
    Ellipses{i} = Ellipses{first - 1};
end

for i=first+roundedHalfWay:last
    Ellipses{i} = Ellipses{last + 1};
end

%%Save Ellipses
save(s,'Ellipses')

end