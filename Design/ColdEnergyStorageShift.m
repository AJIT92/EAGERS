function newDemandC = ColdEnergyStorageShift(DemandE,DemandC, Chiller, TES, grid)
steps = length(DemandE);
Ts = 8760/steps;
TESsize=0;
for i = 1:1:length(TES)
    TESsize = TESsize + TES(i).SysSize;
end
SumOn = find(grid.summerRateTable(2,:)-1,1,'first');
SumOff = find(grid.summerRateTable(2,:)-1,1,'last');
WinOn = find(grid.winterRateTable(2,:)-1,1,'first');
WinOff = find(grid.winterRateTable(2,:)-1,1,'last');
monthDays = [31 28 31 30 31 30 31 31 30 31 30 31];
SummerStart = sum(monthDays(1:grid.summerStartDate(1)-1))+grid.summerStartDate(2);
SummerEnd = sum(monthDays(1:grid.summerEndDate(1)-1))+grid.summerEndDate(2);

ChillerSizesE= zeros(length(Chiller),1);
COPcurve = zeros(length(Chiller(1).COPcurve(:,1)),2,length(Chiller));
for g = 1:1:length(Chiller)
    if strcmp(Chiller(g).ChillType,'electric')
        ChillerSizesE(g) = Chiller(g).SysSize;
        COPcurve(:,1:2,g) = Chiller(g).COPcurve(:,1:2);
    end
end
ChillSize = sum(ChillerSizesE');

%% estimate absorption chiller & modify DemandC
%% Need to add this!!!

%% move peaks with cold energy storage
dayOfWeek = 0;
for day = 1:1:365
    dayOfWeek = dayOfWeek+1;
    if dayOfWeek>7
        dayOfWeek=1;
    end
    if dayOfWeek>=2 && dayOfWeek<=6 && day>=SummerStart && day<SummerEnd && ~isempty(SumOn) % summer peak
        onPeakIndex = (1+(SumOn-1)/Ts+24/Ts*(day-1):(SumOff)/Ts+24/Ts*(day-1))';
        if day ==1
            offPeakIndex = (1:(SumOn-1)/Ts+24/Ts*(day-1))';
        else offPeakIndex = (1+(SumOff)/Ts+24/Ts*(day-2):(SumOn-1)/Ts+24/Ts*(day-1))';
        end
    elseif dayOfWeek>=2 && dayOfWeek<=6 && ~isempty(WinOn) %winter peak
        onPeakIndex = (1+(WinOn-1)/Ts+24/Ts*(day-1):(WinOff)/Ts+24/Ts*(day-1))';
        if day ==1
            offPeakIndex = (1:(WinOn-1)/Ts+24/Ts*(day-1))';
        else offPeakIndex = (1+(WinOff)/Ts+24/Ts*(day-2):(WinOn-1)/Ts+24/Ts*(day-1))';
        end
    else %weekends / no peak pricing
         onPeakIndex = (1+(9-1)/Ts+24/Ts*(day-1):1+(21)/Ts+24/Ts*(day-1))';
        if day ==1
            offPeakIndex = (1:(9-1)/Ts+24/Ts*(day-1))';
        else offPeakIndex = (1+(21)/Ts+24/Ts*(day-2):(9-1)/Ts+24/Ts*(day-1))';
        end
    end
    %Use TES to meet on-peak cooling demand
    OnPeakCool= sum(DemandC(onPeakIndex))*Ts; %total on peak cooling (kWh)
    OffPeakCool= sum(DemandC(offPeakIndex))*Ts;  %total off peak cooling (kWh)
    if (OffPeakCool+min(TESsize,OnPeakCool))/(length(offPeakIndex)*Ts)>ChillSize %ensures that sufficient chilling capcity exists to charge TES in the off-peak period
        TESfill = ((length(offPeakIndex)*Ts)*ChillSize-OffPeakCool)/TESsize; % maximum TES fill based on chiller capacity and off-peak duration
        DemandC(onPeakIndex) = (OnPeakCool-TESfill*TESsize)/(length(onPeakIndex)*Ts);
    elseif OnPeakCool>TESsize
        DemandC(onPeakIndex) = (OnPeakCool-TESsize)/(length(onPeakIndex)*Ts);
        TESfill = 1; %maximum TES fill (sufficient chiller capacity and on-peak demand)
    else DemandC(onPeakIndex) = 0;
        TESfill = OnPeakCool/TESsize; %maximum TES fill (sufficient chiller capacity during off-peak
    end
    %% Flatten the off-peak cooling demand as much as possible
    %%%estimate the electric chiller COP if baseloaded during off-peak
    PowerFracChill = min(1,((TESfill*TESsize+OffPeakCool)/(length(offPeakIndex)*Ts))/ChillSize);
    EstChillCOP = 0;
    for g = 1:1:length(Chiller)
        if strcmp(Chiller(g).ChillType,'electric')
            ChillerCOP = interp1(COPcurve(:,1,g),COPcurve(:,2,g),PowerFracChill)*Chiller(g).COP;
            EstChillCOP = (EstChillCOP*(sum(ChillerSizesE(1:g))-Chiller(g).SysSize)+ChillerCOP*Chiller(g).SysSize)/sum(ChillerSizesE(1:g));
        end
    end
    
    if (TESfill*TESsize+OffPeakCool)>0 %checks that there is a building cooling demand to be met by electric chiller/TES
        
        ChillLoadE = (TESfill*TESsize+OffPeakCool)/EstChillCOP; %convert cooling kWh to kWhe
        avgChillLoadE = ChillLoadE/(length(offPeakIndex)*Ts);%find average power to chillers in kWe
        
        error=TESsize;
        upLimit = mean(DemandE(offPeakIndex))+1.001*avgChillLoadE;
        lowLimit = min(DemandE(offPeakIndex))+avgChillLoadE;
        while error>(.0001*TESsize)
            findThres = linspace(lowLimit,upLimit,10)'*linspace(1,1,length(offPeakIndex)); %kWe threshold to turn on chillers below
            difference = findThres-(DemandE(offPeakIndex)*linspace(1,1,10))';
            areaTest = sum(difference.*(difference>0),2)*Ts; %kWhe under threshold kWe
            if max(areaTest)<ChillLoadE
                upLimit = upLimit +(10/9*ChillLoadE-max(areaTest))/(sum(DemandE(offPeakIndex)<upLimit)*Ts) ; %Raise the initial upper limit slightly for special cases
            else
                I = 1;
                while areaTest(I)<ChillLoadE
                    I = I+1;
                end
                Threshold = findThres(I,1);
                HighDemandE = DemandE(offPeakIndex)>Threshold;
                Threshold = Threshold - (areaTest(I)-ChillLoadE)/(sum(1-HighDemandE)*Ts);
                HighDemandE = DemandE(offPeakIndex)>Threshold;
                error = sum((Threshold-DemandE(offPeakIndex)).*(1-HighDemandE))*Ts-ChillLoadE;
                upLimit = Threshold;
                lowLimit = findThres(max(1,I-1),1);
                if lowLimit>=upLimit && I>=3
                    lowLimit = min(Threshold,findThres((I-2),1));
                end
                if lowLimit>=upLimit && I>=3
                    error = 0;
                end
            end
        end
        %%%elimate dips when building shuts down & still on-peak hours
        if day>1
            A = DemandE(offPeakIndex(1)-10/Ts:offPeakIndex(1)-1);
        if min(A)<Threshold
            j = offPeakIndex(1)-10/Ts-1+find(A,1,'first');
            k = offPeakIndex(length(offPeakIndex));
            while DemandE(j)<Threshold && j > offPeakIndex(1)-5/Ts;
                j = j-1;
            end
            HighDemandE = DemandE(j:k)>Threshold;
            error=sum((Threshold-DemandE(j:k)).*(1-HighDemandE))*Ts-ChillLoadE;
            upLimit = Threshold;
            lowLimit = min(DemandE(j:k))+ChillLoadE/(sum(1-HighDemandE)*Ts);
            while error>(.0001*TESsize)
                findThres = linspace(lowLimit,upLimit,10)'*linspace(1,1,k-j+1); %kWe threshold to turn on chillers below
                difference = findThres-(DemandE(j:k)*linspace(1,1,10))';
                areaTest = sum(difference.*(difference>0),2)*Ts; %kWhe under threshold kWe
                I = 1;
                while areaTest(I)<ChillLoadE
                    I = I+1;
                end
                Threshold = findThres(I);
                HighDemandE = DemandE(j:k)>Threshold;
                Threshold = Threshold - (areaTest(I)-ChillLoadE)/(sum(1-HighDemandE)*Ts);
                HighDemandE = DemandE(j:k)>Threshold;
                error = sum((Threshold-DemandE(j:k)).*(1-HighDemandE))*Ts-ChillLoadE;
                upLimit = Threshold;
                lowLimit = findThres(max(1,I-1),1);
                if lowLimit>=upLimit && I>=3
                    lowLimit = min(Threshold,findThres((I-2),1));
                end
                if lowLimit>=upLimit && I>=3
                    error = 0;
                end
            end
            offPeakIndex = (j:k);
        end
        end
        DemandC(offPeakIndex) = (1-HighDemandE).*(Threshold-DemandE(offPeakIndex))*EstChillCOP;
    end
end
newDemandC = DemandC;
