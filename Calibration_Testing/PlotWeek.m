result = Project.Result.eOut;
StartDate = 1;
Ts = 8760/length(Project.Building(1).DemandE);
month = [0 31 28 31 30 31 30 31 31 30 31 30 31];
monthDays = [0 31 59 90 120 151 181 212 243 273 304 334 365];
monthLabel = ['January  '; 'February '; 'March    '; 'April    '; 'May      '; 'June     '; 'July     '; 'August   '; 'September'; 'October  ';'November ' ;'December ';];
day1 = 1+7*(StartDate-1);
    lastDay = 1+7*StartDate;
    StartMonth = find(monthDays>day1,1) - 1;
    StartDay = day1-monthDays(StartMonth);
    plotAxis = datenum([2014*ones(8,1) StartMonth*ones(8,1)  (StartDay:StartDay+7)'  0*ones(8,1)  0*ones(8,1) 0*ones(8,1)]);
    if StartDay+6<=month(StartMonth+1)
        xlabel(monthLabel(StartMonth,:))
    else 
        xlabel(strcat(monthLabel(StartMonth,:), ' / ', monthLabel(StartMonth+1,:)))
    end
X = datenum(2014,1,1,0,0,0)+((day1-1)+(Ts/24):(Ts/24):lastDay);
PlotType = 'plot';
Y = Project.Building(1).DemandE(1+24/Ts*(day1-1):24/Ts*lastDay);
      

TESshift = result.DemandEshift(1+24/Ts*(day1-1):24/Ts*lastDay);
genPow = result.GenPower(1+24/Ts*(day1-1):24/Ts*lastDay);
        
CO2 = CA.CO2profile;
for i=0:1:(length(CO2)-1)
    CO2new(i*4+1:i*4+4) = CO2(i+1);
end

CO2plot = CO2new((1+24/Ts*(day1-1):24/Ts*lastDay));
[AX,H1,H2] = plotyy(X,TESshift,X,CO2plot);
set(AX(1),'XTickLabel',[{'Mon'},{'Tue'},{'Wed'},{'Thu'},{'Fri'},{'Sat'},{'Sun'}]);
set(AX(2),'XTickLabel',[{'Mon'},{'Tue'},{'Wed'},{'Thu'},{'Fri'},{'Sat'},{'Sun'}]);
set(AX(1),'FontSize',14)
set(AX(2),'FontSize',14)