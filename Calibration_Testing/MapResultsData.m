%% Plot Data on map
%% Map region (1 = lower48, 2 = New England, 3 = South, 4 = Midwest, 5 = West, 6 = Alaska/Hawaii)
MapRegion = 1;

statesAll = {'Alabama';'Alaska';'Arizona';'Arkansas';'California';'Colorado';'Connecticut';'Delaware';'Florida';'Georgia';'Hawaii';'Idaho';'Illinois';'Indiana';'Iowa';'Kansas';
             'Kentucky';'Louisiana';'Maine';'Maryland';'Massachusetts';'Michigan';'Minnesota';'Mississippi';'Missouri';'Montana';'Nebraska';'Nevada';'New Hampshire';'New Jersey';
             'New Mexico';'New York';'North Carolina';'North Dakota';'Ohio';'Oklahoma';'Oregon';'Pennsylvania';'Rhode Island';'South Carolina';'South Dakota';'Tennessee';'Texas';
             'Utah';'Vermont';'Virginia';'Washington';'West Virginia';'Wisconsin';'Wyoming';};

% Y = linspace(1,50,50);
% units = 'Total Demand (GW)';


if MapRegion ==1
%     States2Plot = linspace(1,50,50);
    States2Plot = [1, linspace(3,10,8), linspace(12,50,39)]; %omit hawaii and alaska
elseif MapRegion ==2
    States2Plot = [7, 8, 19, 20, 21, 29, 30, 32, 38, 39, 45];
elseif MapRegion ==3
    States2Plot = [1, 4, 9, 10, 17, 18, 24, 33, 40, 42, 43, 46, 48];    
elseif MapRegion ==4
    States2Plot = [13, 14, 15, 16, 22, 23, 25, 27, 34, 35, 36, 41, 49];    
elseif MapRegion ==5
    States2Plot = [3, 5, 6, 12, 26, 28, 31, 37, 44, 47, 50];  
elseif MapRegion ==6
    States2Plot = [2,11];  
end
states = statesAll(States2Plot);

States2Plot= [];
for i = 1:1:length(states)
    States2Plot(end+1) = find(strcmp(states(i),statesAll));
end
States2Plot = sort(States2Plot);
states = statesAll(States2Plot);
B = sort(Y,'descend');
ColorScaleMax = B(1);
if Y(11) == ColorScaleMax
    ColorScaleMax = B(2);
end
ColorScaleMin = min(0,B(end));
if Y(11) == ColorScaleMin
    ColorScaleMin = B(end-1);
end

if ColorScaleMax<15
    ColorScaleMax = ceil(ColorScaleMax);
elseif ColorScaleMax<=100
    ColorScaleMax = 10*ceil(ColorScaleMax/10);
else ColorScaleMax = 1000*ceil(ColorScaleMax/1000);
end
if ColorScaleMin>-10
    ColorScaleMin = floor(ColorScaleMin);
elseif ColorScaleMin>-40
    ColorScaleMin = 10*floor(ColorScaleMin/10);
elseif ColorScaleMin>=-100
    ColorScaleMin = 100*floor(ColorScaleMin/100);
else ColorScaleMin = 1000*floor(ColorScaleMin/1000);
end
NationalMap(Y(States2Plot),states,units,ColorScaleMax,ColorScaleMin)
