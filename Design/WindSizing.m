function WindSizing(Index,ScaleTo,Repeat,NewSize)
global Project

AnnualPower = zeros(length(Project.Renewable.Wind),1);
for i = 1:1:length(Project.Renewable.Wind)
    wind = Project.Renewable.Wind(i);
    rho = 1.1798-1.3793e-4*wind.Elev+5.667e-9*wind.Elev^2;
    CorrectedWind = (wind.Height/80)^.2*wind.Wind;
    Wind = CorrectedWind.*(wind.Wind>wind.CutIn).*(wind.Wind<wind.ShutDown);
    Wind = min(Wind,wind.Rated);
    AnnualPower(i) = 8760*wind.Eff*.5*rho*(3.1416*(wind.Diam)^2/4)*sum(Wind.^3)/length(Wind)*(1/1000);
end
wind = Project.Renewable.Wind;
if strcmp(ScaleTo,'demand')
    Annual = 8760*sum(Project.Building.DemandE)/length(Project.Building.DemandE);
elseif strcmp(ScaleTo,'generation')
    CHPsize = 0;
    for i = 1:1:length(Project.System.CHP)
        CHPsize = CHPsize + Project.System.CHP(i).SysSize(1);
    end
    Annual = 8760*CHPsize;
end
if strcmp(Repeat,'all')
    oldSize = sum(AnnualPower)/Annual*100;
else oldSize = AnnualPower(Index)/Annual*100;
end
if strcmp(Repeat,'all')
    for i = 1:1:length(wind)
        wind(i).Size = wind(i).Size*(NewSize/oldSize);
        rho = 1.1798-1.3793e-4*wind(i).Elev+5.667e-9*wind(i).Elev^2;
        wind(i).Diam = wind(i).Size/wind(i).Rated.^3*(8/(3.1416*rho*wind(i).Eff));
        wind(i).DemFrac = wind(i).DemFrac*(NewSize/oldSize);
        wind(i).GenFrac = wind(i).GenFrac*(NewSize/oldSize);
    end
else
    wind(Index).Size = wind(Index).Size*(NewSize/oldSize);
    rho = 1.1798-1.3793e-4*wind(Index).Elev+5.667e-9*wind(Index).Elev^2;
    wind(Index).Diam = wind(Index).Size/wind(Index).Rated.^3*(8/(3.1416*rho*wind(Index).Eff));
    wind(Index).DemFrac = wind(Index).DemFrac*(NewSize/oldSize);
    wind(Index).GenFrac = wind(Index).GenFrac*(NewSize/oldSize);
end
Project.Renewable.Wind = wind;