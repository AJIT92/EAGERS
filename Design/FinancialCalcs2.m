function out=FinancialCalcs2(Baseline,Dispatch,grid,Econ,varargin)
% First do the grid analysis for the baseline demand. Scale to the average 
% state electric rate of each month. Then use those scaled electric rates on 
% the dispatched generator grid demand. Then apply the financial calculations.

if ~isempty(varargin)
    state = varargin{1};
    stateName = {'Alabama';'Alaska';'Arizona';'Arkansas';'California';'Colorado';'Connecticut';'Delaware';'Florida';'Georgia';
             'Hawaii';'Idaho';'Illinois';'Indiana';'Iowa';'Kansas';'Kentucky';'Louisiana';'Maine';'Maryland';
             'Massachusetts';'Michigan';'Minnesota';'Mississippi';'Missouri';'Montana';'Nebraska';'Nevada';'NewHampshire';'NewJersey';
             'NewMexico';'NewYork';'NorthCarolina';'NorthDakota';'Ohio';'Oklahoma';'Oregon';'Pennsylvania';'RhodeIsland';'SouthCarolina';
             'SouthDakota';'Tennessee';'Texas';'Utah';'Vermont';'Virginia';'Washington';'WestVirginia';'Wisconsin';'Wyoming';'USavg';};
    stateAbrev = {'AL';'AK';'AZ';'AR';'CA';'CO';'CT';'DE';'FL';'GA';'HI';'ID';'IL';'IN';'IA';'KS';'KY';'LA';'ME';'MD';'MA';'MI';'MN';'MS';'MO';
                  'MT';'NE';'NV';'NH';'NJ';'NM';'NY';'NC';'ND';'OH';'OK';'OR';'PA';'RI';'SC';'SD';'TN';'TX';'UT';'VT';'VA';'WA';'WV';'WI';'WY';'USavg';};  
    I = find(strcmp(state,stateName),1);
    state = char(stateAbrev(I));

    if length(varargin)>1
        user = varargin{2};
    else user = 'Commercial';
    end
    if length(varargin)>2
        natGas = varargin{3};
    end
    if length(varargin)>3
        scale = varargin{4};
    end
else user = 'Commercial';
    state = 'USavg';
    natGas =[];
    scale =1;
end
%% Financial Calculations
startYear = 2014;
Inflation = (1+Econ.Inflation/100);
Years = 2014:1:2033;
inflationPriceFactor=(Inflation.^(Years-startYear));

tmx = [ones(8760,1)*[startYear 1 1]  [0:8759]' ones(8760,1)*[0 0]];
date = datenum(tmx);
month = datevec(date);
month = month(:,2);
[gridRate, sellback]= ElecRate(grid);

%% Demand charges 
DemCharge = zeros(12,1);
GridElec = max(Baseline.Elec,0);
SellBack = max(-Baseline.Elec,0);
for i = 1:12
    [Y, I] = max(GridElec.*(month==i));
    if i>=grid.summerStartDate(1) && i<=grid.summerEndDate(1)
        if gridRate.SumOn(I)>0
            DemCharge(i) = Y*grid.demandCharges(1);
        elseif gridRate.SumMid(I)>0
            DemCharge(i) = Y*grid.demandCharges(2);
        elseif gridRate.SumOff(I)>0
            DemCharge(i) = Y*grid.demandCharges(3);
        end
    else
        if gridRate.WinOn(I)>0
            DemCharge(i) = Y*grid.demandCharges(4);
        elseif gridRate.WinMid(I)>0
            DemCharge(i) = Y*grid.demandCharges(5);
        elseif gridRate.WinOff(I)>0
            DemCharge(i) = Y*grid.demandCharges(6);
        end
    end
end

%% monthly electric
MonthUseCost = zeros(12,1);
for i = 1:12
    MonthUseCost(i) = sum(gridRate.CentskWh.*GridElec.*(month==i)) - sum(sellback.CentskWh.*SellBack.*(month==i));
end

if scale==1%%%Scale monthly bill to state average price
    StateMonthAvg = ProjectUtilityCost('electric',user,Years(end)-startYear+1,state)/100;
    GridPrice = StateMonthAvg;
    MonthRate = zeros(12,1);
    for i = 1:12
        MonthRate(i) = (MonthUseCost(i)+DemCharge(i))/sum(GridElec.*(month==i));
        MonthScaleFactor(i,1) = StateMonthAvg(i)/MonthRate(i);
    end
    for i = 1:1:Years(end)-startYear+1
        gridEscalator(i) = mean(StateMonthAvg(12*(i-1)+1:12*i))/mean(StateMonthAvg(1:12));
    end
else %% use loaded rate structure
    MonthScaleFactor(1:12,1) =1;
    gridEscalator(1) = 1;
    for i = 2:1:Years(end)-startYear+1
        gridEscalator(i) = gridEscalator(i-1)*1.02;
    end
    for i = 1:1:length(Years)
        for j = 1:1:12
            GridPrice(j+12*(i-1)) = (DemCharge(j)+MonthUseCost(j))/sum(Baseline.Elec.*(month==i))*gridEscalator(i);
        end
    end
end
out.Baseline.TotalDemandCharges = sum(DemCharge.*MonthScaleFactor);
out.Baseline.TotalUseCharges = sum(MonthUseCost.*MonthScaleFactor);
out.Baseline.TotalGridBill = out.Baseline.TotalDemandCharges+out.Baseline.TotalUseCharges;
out.GridPriceMonthly = GridPrice;

%% Dispatch Grid Cost
GridElec = max(Dispatch.Elec,0);
SellBack = max(-Dispatch.Elec,0);
for i = 1:12
    [Y, I] = max(GridElec.*(month==i));
    if i>=grid.summerStartDate(1) && i<=grid.summerEndDate(1)
        if gridRate.SumOn(I)>0
            DemCharge(i) = Y*grid.demandCharges(1);
        elseif gridRate.SumMid(I)>0
            DemCharge(i) = Y*grid.demandCharges(2);
        elseif gridRate.SumOff(I)>0
            DemCharge(i) = Y*grid.demandCharges(3);
        end
    else
        if gridRate.WinOn(I)>0
            DemCharge(i) = Y*grid.demandCharges(4);
        elseif gridRate.WinMid(I)>0
            DemCharge(i) = Y*grid.demandCharges(5);
        elseif gridRate.WinOff(I)>0
            DemCharge(i) = Y*grid.demandCharges(6);
        end
    end
end

%%%Scale monthly bill
for i = 1:12
    MonthUseCost(i) = sum(gridRate.CentskWh.*GridElec.*(month==i)) - sum(sellback.CentskWh.*SellBack.*(month==i));
end
out.Dispatch.TotalDemandCharges = sum(DemCharge.*MonthScaleFactor);
out.Dispatch.TotalUseCharges = sum(MonthUseCost.*MonthScaleFactor);
out.Dispatch.TotalGridBill = out.Dispatch.TotalDemandCharges+out.Dispatch.TotalUseCharges;

%% Annual Costs & revenue

%choose either state specific based fuel price or EIA projection
if scale==1%%%Scale monthly bill to state avg
    FuelPrice = ProjectUtilityCost('gas',user,Years(end)-startYear+1,state);
    for i = 1:1:Years(end)-startYear+1
        gasEscalator(i) = mean(FuelPrice(12*(i-1)+1:12*i))/mean(FuelPrice(1:12));
    end
else %% use preloaded gas rates
    AnnualPrice = natGas.fuelPrice(1:length(Years))';
    gasEscalator = (AnnualPrice/AnnualPrice(1));
    for i = 1:1:length(Years)
        FuelPrice(1+12*(i-1):12*i) = AnnualPrice(i);
    end
end
out.FuelPriceMonthly = FuelPrice;
%EIA projection
%    gasEscalator = natGas.fuelPrice';
%    for i = 1:1:Years(end)-startYear+1
%        FuelPrice(12*(i-1)+1:12*i) = gasEscalator(i);
%    end

BaselineFuel = sum(Baseline.Fuel.*(FuelPrice(1:12)/293.01))*gasEscalator;%converts million BTU to kWh
BaselineGridBill = out.Baseline.TotalGridBill*gridEscalator;
BaselineUseCharges = out.Baseline.TotalUseCharges*gridEscalator;
BaselineDemCharges = out.Baseline.TotalDemandCharges*gridEscalator;

FuelCosts = sum(Dispatch.Fuel.*(FuelPrice(1:12)/293.01))*gasEscalator;%converts million BTU to kWh
NewGridBill = out.Dispatch.TotalGridBill*gridEscalator;
NewUseCharges = out.Dispatch.TotalUseCharges*gridEscalator;
NewDemCharges = out.Dispatch.TotalDemandCharges*gridEscalator;

%%%stack replacement
StackReplace = zeros(length(Years),1);
for i = 1:1:length(Years)
    lifespan = Econ.LifekWh/Dispatch.ElecTotProd;
%     if mod(i,Econ.LifeYrs)<1
    if mod(i,lifespan)<1
        StackReplace(i) = Econ.StackReplaceCost*Dispatch.SysSize;
    else StackReplace(i) = 0;
    end
end

%%% Chiller costs
%%Electric Chiller
ChillCost1 = Dispatch.ChillerSize(1)*Econ.ElecChill/3.517; %convert ton to kW

%%absorption chiller, source: Potential use of cold thermal energy storage systems for better efficiency and cost effectiveness, Energy and Buildings, Volume 42 issue 12, december 2010, 2296-2303, M.A. Ehyaei
% cost $/kW= a*ln(size in kW)+b
%small: a = -81.552, b = 778.95
%large: a = -35.469, b = 431.41
if Econ.CurveFit ==1
    %use curve fit
    if Dispatch.ChillerSize(2)>1000
        ChillCost2 = Dispatch.ChillerSize(2)*(-35.469*log(Dispatch.ChillerSize(2))+431.41);
    elseif Dispatch.ChillerSize(2)>0
        ChillCost2 = Dispatch.ChillerSize(2)*(-81.552*log(Dispatch.ChillerSize(2))+778.95);
    else ChillCost2 = 0;
    end
else ChillCost2 = Dispatch.ChillerSize(2)*Econ.AbsorbChill/3.517; %convert ton to kW
end

%%Storage costs
StorageCost = Dispatch.TESsize*Econ.ColdStore/3.517+Dispatch.BatterySize*Econ.Battery; %convert ton to kW
RenewableCost = Dispatch.RenewableSize(1)*Econ.Solar + Dispatch.RenewableSize(2)*Econ.Wind;

NonFCcost = round(StorageCost+ChillCost2+ChillCost1+RenewableCost);
SystemCost = (Econ.InstallCost-Econ.Incentive)*Dispatch.SysSize + NonFCcost;

%%Operation and Maintenance Costs
Total_OM = Econ.CHP_OM*Dispatch.SysSize + Econ.ChillOM*Dispatch.ChillerSize(1) + Econ.AbChillOM*Dispatch.ChillerSize(2) + Econ.ColdOM*Dispatch.TESsize...
    + Econ.BatteryOM*Dispatch.BatterySize + Econ.SolarOM*Dispatch.RenewableSize(1) + Econ.WindOM*Dispatch.RenewableSize(2);
YearlyOandM = inflationPriceFactor*Total_OM+inflationPriceFactor.*StackReplace';

InterestMonth = Econ.Interest/12/100;
months = Econ.FinanceYrs*12;
debtPaymentYearly=(12*(InterestMonth*SystemCost/(1-(1+InterestMonth)^(-months))))*(Years<startYear+Econ.FinanceYrs);

%%% Discounted Payback Period Calculation
%%%Assumes entire installation cost is paid up front
AnnualSavings = (BaselineGridBill+BaselineFuel)-(YearlyOandM+FuelCosts+NewGridBill);
DCF = zeros(length(Years),1);
CumNPV = zeros(length(Years),1);
if sum(AnnualSavings)>0
    zero = [];
    for i = 1:1:length(Years)-1
        DCF(i) = NPV(Inflation,[zero AnnualSavings(i)]);
        CumNPV(i) = sum(DCF(1:i));
        zero = [zero 0];
    end
    if CumNPV(end-1)<SystemCost
        Payback = Years(end)-Years(1);
    else
        i = 1;
        while CumNPV(i)<SystemCost && i+1<length(CumNPV)
            i = i+1;
        end
        if i==1
            Payback = SystemCost/DCF(i);
        else Payback = (i-1)+(SystemCost-CumNPV(i-1))/DCF(i);
        end
    end
else
    Payback = Years(end)-Years(1);
end

out.InstalledCHPcost = SystemCost;
out.CostPerYear = debtPaymentYearly+YearlyOandM+FuelCosts+NewGridBill;
out.CostPerYearkW = out.CostPerYear/Dispatch.SysSize;
out.OandMFinance = debtPaymentYearly+YearlyOandM;
out.NonFCfinanceNPV = NPV(Inflation,(12*(InterestMonth*NonFCcost/(1-(1+InterestMonth)^(-months))))*(Years<startYear+Econ.FinanceYrs));

out.NPVbaselineGridCost=NPV(Inflation,BaselineGridBill);
out.NPVbaselineUseCharges=NPV(Inflation,BaselineUseCharges);
out.NPVbaselineDemCharges=NPV(Inflation,BaselineDemCharges);
out.NPVbaselineFuelCost=NPV(Inflation,BaselineFuel);
out.NPVbaselineOandM = NPV(Inflation,zeros(length(Years),1));
out.NPVbaselineFinance = NPV(Inflation,zeros(length(Years),1));
out.NPVbaselineOandMandFinance = out.NPVbaselineFinance+out.NPVbaselineOandM;


out.NPVnewGridCost=NPV(Inflation,NewGridBill);
out.NPVnewUseCharges=NPV(Inflation,NewUseCharges);
out.NPVnewDemCharges=NPV(Inflation,NewDemCharges);
out.NPVnewFuelCost=NPV(Inflation,FuelCosts);
out.NPVnewOandM = NPV(Inflation,YearlyOandM);
out.NPVnewFinance = NPV(Inflation,debtPaymentYearly);
out.NPVnewOandMandFinance = out.NPVnewOandM+out.NPVnewFinance;
out.NPVnetSavingsFuelandGrid=out.NPVbaselineFuelCost+out.NPVbaselineGridCost-out.NPVnewFuelCost-out.NPVnewGridCost;

out.actualYears = Years;
out.BaselineGridBill = BaselineGridBill;
out.BaselineFuel = BaselineFuel;
out.FuelCosts = FuelCosts;
out.NewGridBill = NewGridBill;
out.Payback = Payback;

out.Year1baseCharges = [out.Baseline.TotalDemandCharges out.Baseline.TotalUseCharges BaselineFuel(1) ];
out.Year1dispatchCharges = [out.Dispatch.TotalDemandCharges out.Dispatch.TotalUseCharges FuelCosts(1)];

function value=NPV(irr,cashflows)

value=0;
for i=1:length(cashflows)
    value=value+cashflows(i)/(irr)^i;
end