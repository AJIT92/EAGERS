function SparkSpreadMap(user,years,color)
global MapHandle
MapHandle =figure(1);
listUser = {'Residential'; 'Commercial';'Industrial';'Electric Utility';};

statesAll = {'Alabama';'Alaska';'Arizona';'Arkansas';'California';'Colorado';'Connecticut';'Delaware';'Florida';'Georgia';
             'Hawaii';'Idaho';'Illinois';'Indiana';'Iowa';'Kansas';'Kentucky';'Louisiana';'Maine';'Maryland';
             'Massachusetts';'Michigan';'Minnesota';'Mississippi';'Missouri';'Montana';'Nebraska';'Nevada';'New Hampshire';'New Jersey';
             'New Mexico';'New York';'North Carolina';'North Dakota';'Ohio';'Oklahoma';'Oregon';'Pennsylvania';'Rhode Island';'South Carolina';
             'South Dakota';'Tennessee';'Texas';'Utah';'Vermont';'Virginia';'Washington';'West Virginia';'Wisconsin';'Wyoming';};
         
% load(fullfile(Model_dir,'component library', 'NatGas','RateData','GasRate'))
% load(fullfile(Model_dir,'component library', 'Grid','RateData','ElecRate'))
user2 = listUser{user};
stateAbrev = {'AL';'AK';'AZ';'AR';'CA';'CO';'CT';'DE';'FL';'GA';'HI';'ID';'IL';'IN';'IA';'KS';'KY';'LA';'ME';'MD';'MA';'MI';'MN';'MS';'MO';
              'MT';'NE';'NV';'NH';'NJ';'NM';'NY';'NC';'ND';'OH';'OK';'OR';'PA';'RI';'SC';'SD';'TN';'TX';'UT';'VT';'VA';'WA';'WV';'WI';'WY';};
GasProjection = zeros(50,12*years);
ElecProjection = zeros(50,12*years);
for s = 1:1:50
    state = char(stateAbrev(s));
    GasProjection(s,:) = ProjectUtilityCost('gas',user2,years,state);
    ElecProjection(s,:) = ProjectUtilityCost('electric',user,years,state);
end
SparkSpread = sum(ElecProjection/100-GasProjection/293.297,2)/(12*years);%average spark spread in each state over the next ___years
NatSparkSpread = sum(ProjectUtilityCost('electric',user2,years,'USavg')/100-ProjectUtilityCost('gas',user2,years,'USavg')/293.297,2)/(12*years);
RangeUpper = max([SparkSpread(1); SparkSpread(3:10);SparkSpread(12:50)])-NatSparkSpread;
RangeLower = NatSparkSpread-min([SparkSpread(1);SparkSpread(3:10);SparkSpread(12:50)]);
for s = 1:1:50
    if SparkSpread(s)>NatSparkSpread
        NormSparkSpread(s) = (SparkSpread(s)-NatSparkSpread)/RangeUpper;
    else
       NormSparkSpread(s) = -(NatSparkSpread-SparkSpread(s))/RangeLower; 
    end
end

units='';
States2Plot= [];
states = statesAll;
A = NormSparkSpread;
for i = 1:1:50
    if strcmp(statesAll(2),states(i)) || strcmp(statesAll(11),states(i))
        if length(states)<=2 %omit hawaii and alaska
            States2Plot(end+1) = find(strcmp(states(i),statesAll));
        end
    else
        States2Plot(end+1) = find(strcmp(states(i),statesAll));
    end
end
States2Plot = sort(States2Plot);
states = statesAll(States2Plot);

[units,ColorScaleMax,ColorScaleMin,ticks] = ScaleAndUnits('',NormSparkSpread(States2Plot));
NationalMapPlot(NormSparkSpread(States2Plot),states,units,ColorScaleMax,ColorScaleMin,ticks,color)

MapHandle =figure(2);
ElecCostNow = sum(ElecProjection(:,1:12),2)/12;
[units,ColorScaleMax,ColorScaleMin,ticks] = ScaleAndUnits('Average Commercial Electric Rate (¢/kWh)',ElecCostNow(States2Plot));
NationalMapPlot(ElecCostNow(States2Plot),states,units,ColorScaleMax,ColorScaleMin,ticks,color)
MapHandle =figure(3);
GasCostNow = sum(GasProjection(:,1:12),2)/12;
[units,ColorScaleMax,ColorScaleMin,ticks] = ScaleAndUnits('Average Commercial Gas Rate ($/mmBTU)',GasCostNow(States2Plot));
NationalMapPlot(GasCostNow(States2Plot),states,units,ColorScaleMax,ColorScaleMin,ticks,color)

function [units,ColorScaleMax,ColorScaleMin,ticks] = ScaleAndUnits(label,A)
units = label;
B = sort(A,'descend');
ColorScaleMax = ceil(B(1));
if A(11) == ColorScaleMax
    ColorScaleMax = ceil(B(2));
end
ColorScaleMin = floor(B(end));
if ColorScaleMax<2
    ColorScaleMax = ceil(ColorScaleMax);
    ticks = 5*ceil(ColorScaleMax-ColorScaleMin)+1;
elseif ColorScaleMax<20
    ColorScaleMax = ceil(ColorScaleMax);
    ticks = ceil(ColorScaleMax-ColorScaleMin)+1;
elseif ColorScaleMax<100
    ColorScaleMax = 10*ceil(ColorScaleMax/10);
    ticks = ceil((ColorScaleMax-ColorScaleMin)/100);
else ColorScaleMax = 1000*ceil(ColorScaleMax/1000);
    ticks = ceil((ColorScaleMax-ColorScaleMin)/1000);
end
if ColorScaleMin>-10
    ColorScaleMin = floor(ColorScaleMin);
elseif ColorScaleMin>-100
    ColorScaleMin = 10*floor(ColorScaleMin/10);
else ColorScaleMin = 1000*floor(ColorScaleMin/1000);
end
if ticks<3
    ticks = ceil(ColorScaleMax-ColorScaleMin)/10+1;
end
if    ticks>10
    ticks = 10;
end