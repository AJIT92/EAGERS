function [Rate, Sellback] = ElecRate(grid)
SummerRate = zeros(8760,1);  
WinterRate = zeros(8760,1);  
Year = 2014;
tmx = [ones(8760,1)*[Year 1 1]  [0:8759]' ones(8760,1)*[0 0]];
date = datenum(tmx);

SummerStart = datenum([num2str(grid.summerStartDate(1)),'/',num2str(grid.summerStartDate(2)),'/',num2str(Year)]);
SummerEnd = datenum([num2str(grid.summerEndDate(1)),'/',num2str(grid.summerEndDate(2)),'/',num2str(Year)]);
summer = date>=SummerStart & date<=(SummerEnd+.99);

Rate1 = (grid.summerRateTable==1);
Rate2 = (grid.summerRateTable==2);
Rate3 = (grid.summerRateTable==3);
SumRates = Rate1*grid.summerPrice(1)+Rate2*grid.summerPrice(2)+Rate3*grid.summerPrice(3);
Rate4 = (grid.winterRateTable==1);
Rate5 = (grid.winterRateTable==2);
Rate6 = (grid.winterRateTable==3);
WinRates = Rate4*grid.winterPrice(1)+Rate5*grid.winterPrice(2)+Rate6*grid.winterPrice(3);

day = weekday(date(1));
for i = 1:1:365
    SummerRate(1+24*(i-1):24*i) = SumRates(day,1:24)';
    WinterRate(1+24*(i-1):24*i) = WinRates(day,1:24)';
    Rate.SumOn(1+24*(i-1):24*i) =  Rate1(day,1:24)';
    Rate.SumMid(1+24*(i-1):24*i) =  Rate2(day,1:24)';
    Rate.SumOff(1+24*(i-1):24*i) =  Rate3(day,1:24)';
    Rate.WinOn(1+24*(i-1):24*i) =  Rate4(day,1:24)';
    Rate.WinMid(1+24*(i-1):24*i) =  Rate5(day,1:24)';
    Rate.WinOff(1+24*(i-1):24*i) =  Rate6(day,1:24)';
    day = day+1;
    if day ==8
        day =1;
    end
end
Rate.CentskWh = summer.*SummerRate + (1-summer).*WinterRate;
if ~isfield(grid,'SellBackPerc');
    grid.SellBackPerc =100;
end
if grid.SellBackRate ==-1
    Sellback.CentskWh = Rate.CentskWh*grid.SellBackPerc/100;
elseif grid.SellBackRate ==0
    Sellback.CentskWh = zeros(8760,1);
else Sellback.CentskWh = (zeros(8760,1)+1)*grid.SellBackRate;
end