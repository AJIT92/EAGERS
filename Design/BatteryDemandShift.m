function newDemandE = BatteryDemandShift(DemandE, Battery)
steps = length(DemandE);
Ts = 8760/steps;

%% Battery Characteristics
Size = zeros(length(Battery),1);
OptSize = zeros(length(Battery),1);
PeakDisch = zeros(length(Battery),1);
PeakCharge = zeros(length(Battery),1);
DischResist = zeros(length(Battery),1);
ChargeResist = zeros(length(Battery),1);
MaxDOD = zeros(length(Battery),1);
Voltage = zeros(length(Battery),1);
VoltCurveX = zeros(length(Battery(1).VoltCurve(:,1)),length(Battery));
VoltCurveY = zeros(length(Battery(1).VoltCurve(:,1)),length(Battery));
for i = 1:1:length(Battery)
    Size(i)         = Battery(i).Size;
    OptSize(i)      = Battery(i).OptSize/100;
    PeakDisch(i)    = Battery(i).PeakDisch;
    PeakCharge(i)   = Battery(i).PeakCharge;
    DischResist(i)  = Battery(i).DischResist;
    ChargeResist(i) = Battery(i).ChargeResist;
    MaxDOD(i)       = Battery(i).MaxDOD;
    Voltage(i)      = Battery(i).Voltage;
    VoltCurveX(:,i) = Battery(i).VoltCurve(:,1);
    VoltCurveY(:,i) = Battery(i).VoltCurve(:,2);
end
DischCurrent = PeakDisch.*Size./Voltage*1000;
DischResistScaled = (100./DischCurrent).*DischResist*(1/1000); %Scale so the loss of power is equivelant to that specified at 100Amps
DischVoltLoss = DischCurrent.*DischResistScaled;
UseableSize  = Size.*(MaxDOD/100).*(Voltage-DischVoltLoss)./Voltage;
totalUseSize = sum(UseableSize);
PeakDischPower = sum(DischCurrent.*(Voltage-DischVoltLoss))/1000;
ChargeCurrent = PeakCharge.*Size./Voltage*1000;
ChargeResistScaled = (100./ChargeCurrent).*ChargeResist*(1/1000); %Scale so the loss of power is equivelant to that specified at 100Amps
ChargeVoltLoss = ChargeCurrent.*ChargeResistScaled;
PeakChargePower = sum(ChargeCurrent.*(Voltage+ChargeVoltLoss))/1000;
ChargeEff = (Voltage-DischVoltLoss)./(Voltage+ChargeVoltLoss);


%% shave peaks with a battery
monthDays = [0 31 59 90 120 151 181 212 243 273 304 334 365];
for month = 1:1:12
    %% find threshold where battery will come on
    monthProf = DemandE(1+monthDays(month)*24/Ts:monthDays(month+1)*24/Ts);
    Threshold = max(monthProf)- PeakDischPower;
    i = max(find(monthProf>Threshold,1),2);
    event = 0;
    MaxBat(1) = 1;
    while i <=length(monthProf)
        j= i+monthDays(month)*24/Ts;
        if (DemandE(j)>Threshold && DemandE(j-1)<=Threshold) || event ==0
            Energy = 0; %kWh of electric demand
            ThreshCross = j;
            event = event+1;
            ChargeCap(event) = 0;
        end
        Energy = Energy+max((DemandE(j)-Threshold)*Ts,0);%%add up energy being shifted
        if Energy>totalUseSize % ensure batterry is big enough
            Threshold = Threshold+(Energy-totalUseSize)/((j-ThreshCross)*Ts);
            Energy = sum(max(DemandE(ThreshCross:j)-Threshold,0)*Ts);
        end
        MinBat(event) = MaxBat(event)-Energy/sum(Size);
        
        if event>0 && DemandE(j)<=Threshold && DemandE(j-1)>Threshold
            next = j;
            while DemandE(next)<Threshold && next<8760/Ts % ensure there will be time to charge the battery
                next = next+1;
            end
            ChargeCap(event) = sum(Threshold -DemandE(j:next-1))*Ts;
            while MinBat(event)<(1-max(MaxDOD/100))
                NewThreshold = Threshold+.25*(max(monthProf)-Threshold);
                NewEnergy = sum(max(DemandE(ThreshCross:j-1)-NewThreshold,0)*Ts);
                if MaxBat(event)-NewEnergy/sum(Size) > 1-max(MaxDOD/100)
                    MinBat(event) = MaxBat(event)-NewEnergy/sum(Size);
                end
                Threshold = NewThreshold;
                ChargeCap(event) = sum(Threshold -DemandE(j:next-1))*Ts;
            end
            MaxBat(event+1) = min(MinBat(event)+ChargeCap(event)/(totalUseSize*sum(ChargeEff.*(UseableSize/totalUseSize))),1);
        end
        i = i+1;
    end
    %% now dispatch the battery to meet all times > threshold
    start = 1+monthDays(month)*24/Ts;
    finish = monthDays(month+1)*24/Ts;
    BatSOC(1:length(Battery))       = 1;
    Volts(1:length(Battery))        = 0;
    CurrentEst(1:length(Battery))   = 0;
    TerminalVolt(1:length(Battery)) = 0;
    Current(1:length(Battery))      = 0;
    t = start;
    while t<finish
        if DemandE(t)>Threshold && min(BatSOC'>(1-MaxDOD./100))>0%use battery
            PowBat = DemandE(t)-Threshold;
            DemandE(t) = Threshold;
            for i = 1:1:length(Battery)
                Volts(i) = interp1(VoltCurveX(:,i),VoltCurveY(:,i),BatSOC(i)*100);
                CurrentEst(i) = PowBat*(UseableSize(i)/totalUseSize)/Volts(i)*1000;
                TerminalVolt(i) = Volts(i)-CurrentEst(i)*DischResistScaled(i);
                Current(i) = PowBat*(UseableSize(i)/totalUseSize)/TerminalVolt(i)*1000;
                BatSOC(i) = BatSOC(i)-(Current(i)*Volts(i)/1000)*Ts/Size(i);
            end
        elseif min(BatSOC)<1%charge battery if capacity exists
            next = t+1;
            while DemandE(next)<=Threshold && next<8760/Ts
                next = next+1;
            end
            Eshift = sum((1-BatSOC)'.*Size);
            MinChargePower = Eshift./(ChargeEff*(next-t)*Ts);
            AvgLoad = mean(DemandE(t:next))+sum(MinChargePower);
            if min(Threshold>DemandE(t:next-1)+sum(MinChargePower))==1 
                ChargePow = sum(MinChargePower);
            else ChargePow = max(min(Threshold,AvgLoad)-DemandE(t),0);
            end
            for i = 1:1:length(Battery)
                Volts(i) = interp1(VoltCurveX(:,i),VoltCurveY(:,i),BatSOC(i)*100);
                if strcmp(Battery(i).ChargeMeth,'smooth') %Smoothest demand profile
                    CurrentEst(i) = ChargePow*(UseableSize(i)/totalUseSize)/Volts(i)*1000;
                    TerminalVolt(i) = Volts(i)+CurrentEst(i)*ChargeResistScaled(i);
                    Current(i) = ChargePow*(UseableSize(i)/totalUseSize)/TerminalVolt(i)*1000;
                elseif strcmp(Battery(i).ChargeMeth,'CC') || strcmp(Battery(i).ChargeMeth,'CV')
                    if strcmp(Battery(i).ChargeMeth,'CC') %Constant current
                        if sum(ChargeCurrent(1:i).*Volts(1:i)/1000)<Threshold-DemandE(t)
                            Current(i) = ChargeCurrent(i);
                        else Current(i) = 0;
                        end
                        TerminalVolt(i) = Volts(i)+Current(i)*ChargeResistScaled(i);
                    elseif strcmp(Battery(i).ChargeMeth,'CV') %Constant voltage
                        TerminalVolt(i) = interp1(VoltCurveX(:,i),VoltCurveY(:,i),100)+.2;% a voltage slightly higher than peak voltage to ensure 100% charging
                        Current(i) = (TerminalVolt(i)-Volts(i))/ChargeResistScaled(i);
                        if sum(Current(1:i).*TerminalVolt(1:i)/1000)>Threshold-DemandE(t)
                            Current(i) = 0;
                        end
                    end
                end
                ChargeFrac(i) = min((1-BatSOC(i))/((Current(i)*Volts(i)/1000)*Ts/Size(i)),1); %prevents overcharging battery
                BatSOC(i) = BatSOC(i)+(Current(i)*Volts(i)/1000)*Ts/Size(i)*ChargeFrac(i);
            end
            PowBat = sum(TerminalVolt.*Current.*ChargeFrac)/1000;
            DemandE(t) = DemandE(t)+PowBat;
        end
        t = t+1;
    end
end
newDemandE = DemandE;
