function autoSizing(param,component,varargin)
global Project

nSys = length(Project.System.CHP);
SysSize = zeros(nSys,1);
for g= 1:nSys
    SysSize(g) = Project.System.CHP(g).SysSize(1,1);
end
CHPsize=sum(SysSize); %Original system size
switch component
    case 'CHP'
        nSys = length(Project.System.CHP);
        TurnDown = zeros(nSys,1);
        for g= 1:nSys
            TurnDown(g) = Project.System.CHP(g).TurnDown;
        end
        CurrentSize = (CHPsize);
        highGuess = max(Project.Building.DemandE);
        lowGuess =  min(Project.Building.DemandE);
    case 'Battery'
        nSys = length(Project.System.Battery);
        Size = zeros(nSys,1);
        for g= 1:nSys
            BatSize(g) = Project.System.Battery(g).Size;
        end
        BatSizeT=sum(BatSize); %Original system size
        CurrentSize = (BatSizeT);
        OrderedDemand = sort(Project.Building.DemandE);
        Guess1 = varargin{1}*(OrderedDemand(floor(.98*length(OrderedDemand)))-CHPsize);% Capacity necessary to meet x hours of the difference between peak generation and peak demand 
        Guess2 = varargin{1}*0.2*(OrderedDemand(floor(.98*length(OrderedDemand))) - OrderedDemand(floor(.02*length(OrderedDemand)))); % Capacity to meet x hours of 20% of the variability in building load
        Guess3 = varargin{1}*CHPsize;
        highGuess = min(Guess3,max(Guess1,Guess2));
        lowGuess =  0;
    case 'TES'
        nSys = length(Project.System.CHP);
        for g= 1:nSys
            TESsize(g) = Project.System.TES(g).SysSize;
        end
        TESsizeT = sum(TESsize);
        CurrentSize = (TESsizeT);
        highGuess = 100;
        lowGuess =  0;
end


if ~strcmp(param,'cost')
    stateName = {'Alabama';'Alaska';'Arizona';'Arkansas';'California';'Colorado';'Connecticut';'Delaware';'Florida';'Georgia';
             'Hawaii';'Idaho';'Illinois';'Indiana';'Iowa';'Kansas';'Kentucky';'Louisiana';'Maine';'Maryland';
             'Massachusetts';'Michigan';'Minnesota';'Mississippi';'Missouri';'Montana';'Nebraska';'Nevada';'NewHampshire';'NewJersey';
             'NewMexico';'NewYork';'NorthCarolina';'NorthDakota';'Ohio';'Oklahoma';'Oregon';'Pennsylvania';'RhodeIsland';'SouthCarolina';
             'SouthDakota';'Tennessee';'Texas';'Utah';'Vermont';'Virginia';'Washington';'WestVirginia';'Wisconsin';'Wyoming';};
    StateNum = find(strcmp(Project.State,stateName),1,'first');
    stateAbrev = {'AL';'AK';'AZ';'AR';'CA';'CO';'CT';'DE';'FL';'GA';'HI';'ID';'IL';'IN';'IA';'KS';'KY';'LA';'ME';'MD';'MA';'MI';'MN';'MS';'MO';
                  'MT';'NE';'NV';'NH';'NJ';'NM';'NY';'NC';'ND';'OH';'OK';'OR';'PA';'RI';'SC';'SD';'TN';'TX';'UT';'VT';'VA';'WA';'WV';'WI';'WY';};
    State = stateAbrev(StateNum);
    [CO2, NOx, SO2] = EmissionProfile(State,0,2);
end
error = CurrentSize;
while error>.005*CurrentSize
    guessSizes = linspace(lowGuess,highGuess,10);
    CostFunction = zeros(length(guessSizes),1);
    i = 1;
    while i<=length(guessSizes)
        switch component
            case 'CHP'
                %scale everything else proportionately
                sol = 0;
                win = 0;
                if isfield(Project,'Renewable')
                    if isfield(Project.Renewable,'Solar')
                        for l = 1:1:length(Project.Renewable.Solar)
                            sol = sol + sum(Project.Renewable.Solar(l).DemFrac);
                        end
                    end
                    if isfield(Project.Renewable,'Wind')
                        for l = 1:1:length(Project.Renewable.Wind)
                            win = win + sum(Project.Renewable.Wind(l).DemFrac);
                        end
                    end
                end
                if isfield(Project.System,'Battery')
                    BatSizeMinutes = Project.System.Battery(1).SizeMin; %Original system size
                else BatSizeMinutes = 0;
                end
                if isfield(Project.System,'Chiller')
                    ChillOptSize = Project.System.Chiller(1).OptSize;
                else ChillOptSize = 0;
                end
                if isfield(Project.System,'TES')
                    TESoptSize = Project.System.TES(1).OptSize;
                else TESoptSize = 0;
                end
                %% scale CHP system
                if strcmp(char(Project.Control.Name),'1_BaseLoad') && Project.Utilities.Grid.SellBackRate ==0
                    ScaleSystemComponents(100,ChillOptSize,TESoptSize,BatSizeMinutes,sol,2,win,2,2)
                    guessSizes = 0;
                    for chp =1:1:length(Project.System.CHP)
                        guessSizes = guessSizes+Project.System.CHP(chp).SysSize(1);
                    end
                else
                    for chp =1:1:length(Project.System.CHP)
                        Project.System.CHP(chp).SysSize(1) = guessSizes(i)*SysSize(chp)/CHPsize;
                        Project.System.CHP(chp).SysSize(2) = guessSizes(i)*SysSize(chp)/CHPsize/TurnDown(chp);
                        OptimalFCsizing(chp,'editMaxPower');
                    end
                    ScaleSystemComponents(Project.System.CHP(1).OptSize,ChillOptSize,TESoptSize,BatSizeMinutes,sol,2,win,2,2)
                end
            case 'Battery'
                for bat =1:1:length(Project.System.Battery)
                    Project.System.Battery(bat).SysSize(1) = guessSizes(i)*BatSize(bat)/BatSizeT;
                end
            case 'TES'
                for tes =1:1:length(Project.System.TES)
                    Project.System.TES(tes).OptSize = guessSizes(i);
                end
                OptimalTESsizing('OptSize')
        end
        RESULT = runAnalyses(Project);
        switch param
            case 'cost'
                CostFunction(i) = (RESULT.costOut.NPVnewUseCharges + RESULT.costOut.NPVnewDemCharges + RESULT.costOut.NPVnewFuelCost +...
                            RESULT.costOut.NPVnewOandM + RESULT.costOut.NPVnewFinance);
            case 'CO2'
                BoilerCO2 = 0.3988; %lb CO2 per kWh
                CostFunction(i) = sum([sum(RESULT.eOut.Grid_purchases_hourly_kWh.*CO2) , sum(RESULT.eOut.CO2) , sum(BoilerCO2*RESULT.eOut.BoilerHeatHour)]);
            case 'SO2'
                BoilerSO2 = 2.5e-3; %lb SO2/kWh
                CostFunction(i) = sum([sum(RESULT.eOut.Grid_purchases_hourly_kWh.*SO2) , sum(RESULT.eOut.SO2) , sum(BoilerSO2*RESULT.eOut.BoilerHeatHour)]);
            case 'NOx'
                BoilerNOx = 7.7e-4; %lbNOx/kWh
                CostFunction(i) = sum([sum(RESULT.eOut.Grid_purchases_hourly_kWh.*NOx) , sum(RESULT.eOut.NOx) , sum(BoilerNOx*RESULT.eOut.BoilerHeatHour)]);
        end
        i = i+1;
    end
    if length(guessSizes)>1
        [minCost, Index] = min(CostFunction);
        lowGuess = guessSizes(max(1,Index-1));
        highGuess = guessSizes(min(length(guessSizes),Index+1));
        if Index ==1 || Index ==length(CostFunction)
            error = 0.01*error;
        else error = abs(highGuess-lowGuess);
        end
    else error = 0;
        Index = 1;
    end
end
newSize = guessSizes(Index);
switch component
    case 'CHP'
        for chp =1:1:length(Project.System.CHP)
            Project.System.CHP(chp).SysSize(1) = newSize*SysSize(chp)/CHPsize;
            Project.System.CHP(chp).SysSize(2) = newSize*SysSize(chp)/CHPsize/TurnDown(chp);
        end
    case 'Battery'
        for bat =1:1:length(Project.System.Battery)
            Project.System.Battery(bat).SysSize(1) = newSize*BatSize(bat)/BatSizeT;
        end
    case 'TES'
        for tes =1:1:length(Project.System.TES)
            Project.System.TES(tes).OptSize = newSize*TESsize(tes)/TESsizeT;
        end
        OptimalTESsizing('OptSize')
end