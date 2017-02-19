function Projection = ProjectUtilityCost(utility,user,years,state)
global Model_dir
if strcmp(utility,'gas')
    load(fullfile(Model_dir,'System Library', 'NatGas','RateData','GasRate'))
    Date1 = GasRate.Date(:,1);
    CityGate = GasRate.(state).CityGate;
    Residential = GasRate.(state).Residential;
    Commercial = GasRate.(state).Commercial;
    Industrial = GasRate.(state).Industrial;
    ForElectricGeneration = GasRate.(state).ForElectricGeneration;

    monthAvg = zeros(1,12);
    for j = 1:1:floor(length(Date1)/12)
        monthAvg(1:12) = monthAvg+CityGate(12*(j-1)+1:12*j);
        annualAvg(j) = mean(CityGate(12*(j-1)+1:12*j));
    end
    monthAvg = monthAvg/j;
elseif strcmp(utility,'electric')
    load(fullfile(Model_dir,'System Library', 'Grid','RateData','ElecRate'))
    Date = ElecRate.Date;
    Residential = ElecRate.(state).Residential;
    Commercial = ElecRate.(state).Commercial;
    Industrial = ElecRate.(state).Industrial;

    monthAvg = zeros(1,12);
    for j = 1:1:floor(length(Date)/12)
        monthAvg(1:12) = monthAvg+Commercial(12*(j-1)+1:12*j);
        annualAvg(j) = mean(Commercial(12*(j-1)+1:12*j));
    end
    monthAvg = monthAvg/j;
    if strcmp(state,'CA')
        annualAvg(1:3) = annualAvg(4); %%accounts for california energy crisis (Enron)
    end
end

X = linspace(j+1,j+years,years);
A = polyfit(linspace(1,j,j),annualAvg,2);
Y = A(1)*X.^2+A(2)*X+A(3);
if A(2)<0 || A(1)<0
    A = polyfit(linspace(1,j,j),annualAvg,1);
    Y = A(1)*X+A(2);
    dif = Y(1)-(annualAvg(end)+A(1));
    Y = Y-dif+dif*sqrt(linspace(.001,1,years));
end

for j = 1:1:length(Y)
    Projection(12*(j-1)+1:12*j) = Y(j)+monthAvg-mean(monthAvg);
end

switch user
    case 'Commercial'
        Projection = Projection - (Projection(1)-Commercial(end));
    case 'Residential'
        Projection = Projection - (Projection(1)-Residential(end));
    case 'Industrial'
        Projection = Projection - (Projection(1)-Industrial(end));
    case 'Electric Utility'
        Projection = Projection - (Projection(1)-ForElectricGeneration(end));
end