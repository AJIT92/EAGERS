function scaleCost = updateGeneratorCost(Time)
global Plant DateSim %DateSim: Current time in the simulation.
nS = length(Time);
Timestamp = Time/24 + DateSim;
Source = {zeros(length(Plant.Generator),1)};
for i = 1:1:length(Plant.Generator)
    if strcmp(Plant.Generator(i).Type,'Utility')
        Source(i) = {Plant.Generator(i).Source};
        utility = Plant.Generator(i).VariableStruct;
        if strcmp(Source(i), 'Electricity')
            day = weekday(Timestamp);
            [Y,M,D] = datevec(Timestamp(1));
            WinStart = datenum([Y,utility.WinStartMonth,utility.WinStartDay]);
            if M>utility.WinStartMonth || (M==utility.WinStartMonth && D>=utility.WinStartDay)
                SumStart = datenum([Y+1,utility.SumStartMonth,utility.SumStartDay]);
            else
                SumStart = datenum([Y,utility.SumStartMonth,utility.SumStartDay]);
            end
            Rate = zeros(nS,1);
            [~,~,~,H,~,~] = datevec(Timestamp);
            H = H+1; %zeroth hour is column 1 in rate matrix
            for t = 1:1:length(H) 
                if Timestamp(t)>=SumStart && Timestamp(t)<WinStart
                    Rate(t) = utility.SumRates(utility.SumRateTable(day(t),H(t),1));
                else
                    Rate(t) = utility.WinRates(utility.WinRateTable(day(t),H(t),1));
                end
            end
            Utility(i).Rate = Rate;
        elseif strcmp(Source(i), 'NG')
            %gas utility set up as daily prices for a leap year
            date1 = datevec(utility.Timestamp(1));
            year1 = date1(1);
            datenow = datevec(Timestamp);
            interpDate = datenum([year1*ones(length(Time),1), datenow(:,2:end)]);
            Utility(i).Rate = interp1(utility.Timestamp,utility.Rate,interpDate)/293.1; %interpolate & convert gas rate from $/MMBTU to $/kWh;
        else
            %% need to add something for district heating/cooling
        end
    end
end

scaleCost = zeros(nS,length(Plant.Generator));
for i = 1:1:length(Plant.Generator)
    if strcmp(Plant.Generator(i).Type,'Utility')
        scaleCost(:,i) = Utility(i).Rate;
    elseif ~isfield(Plant.Generator(i).OpMatB,'Stor') && ~strcmp(Plant.Generator(i).Source, 'Renewable')
        Uindex = find(strcmp(Source,Plant.Generator(i).Source),1);
        if isempty(Uindex)
           scaleCost(:,i) = 1; %no utility (don't scale costs)
        else
            scaleCost(:,i) = Utility(Uindex).Rate;
        end
    end
end