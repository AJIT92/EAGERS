function OptimalBatterySizing(sizeChange)
% determines OptSize from SysSiz & vice versa for different control cases
global Project
demand=Project.Building; %the building load structure
grid=Project.Utilities.Grid;
Battery = Project.System.Battery;
%temporarliy remove battery from system
Project.System = rmfield(Project.System,'Battery');
systems=Project.System;
if isfield(Project,'Renewable')
    renew = Project.Renewable;
else renew = 1;
end
[DemandE, DemandH,out] = nonCHPdispatch(systems,demand,grid,renew); 
steps = length(DemandE);
Ts = 8760/steps;
%Find Max peak that can be reduced
dayProf= zeros(365,1);
for day = 1:1:365
    dayProf(day,1:24/Ts) = DemandE(1+24/Ts*(day-1):24/Ts*day);
end
[HighPeak, Day] = max(max(dayProf,[],2));
PeakMean = mean(dayProf(Day,:));
[HighPeak, tPeak] = max(dayProf(Day,:));
tStart = tPeak;
while tStart>1 && dayProf(Day,tStart-1)>PeakMean
    tStart = tStart-1;
end
tEnd = tPeak;
while tEnd<24/Ts && dayProf(Day,tEnd+1)>PeakMean
    tEnd = tEnd+1;
end
OptimalBatSize = sum(dayProf(Day,tStart:tEnd))-PeakMean*(tEnd-tStart)*Ts;
CHPsize = 0;
for i = 1:1:length(Project.System.CHP)
    CHPsize = CHPsize + Project.System.CHP(i).SysSize(1);
end

Size = zeros(length(Battery),1);
SizeMin = zeros(length(Battery),1);
OptSize = zeros(length(Battery),1);
PeakDisch = zeros(length(Battery),1);
DischResist = zeros(length(Battery),1);
MaxDOD = zeros(length(Battery),1);
Voltage = zeros(length(Battery),1);
for i = 1:1:length(Battery)
    Size(i)         = Battery(i).Size;
    SizeMin(i)      = Battery(i).SizeMin;
    OptSize(i)      = Battery(i).OptSize/100;
    PeakDisch(i)    = Battery(i).PeakDisch;
    DischResist(i)  = Battery(i).DischResist;
    MaxDOD(i)       = Battery(i).MaxDOD;
    Voltage(i)      = Battery(i).Voltage;
end
totalSize = sum(Size);
totalSizeMin = sum(SizeMin/60*CHPsize);
DischCurrent = PeakDisch.*Size./Voltage*1000;
DischResistScaled = (100./max(DischCurrent,.01)).*DischResist*(1/1000); %Scale so the loss of power is equivelant to that specified at 100Amps
DischVoltLoss = DischCurrent.*DischResistScaled;
UseableSize  = Size.*(MaxDOD/100).*(Voltage-DischVoltLoss)./Voltage;
% UseableSize  = (1/PeakDisch).*DischCurrent.*(Voltage-DischVoltLoss).*(MaxDOD/100);
totalUseSize = sum(UseableSize);

for i = 1:1:length(Battery)
    if strcmp(sizeChange,'editSize')
        SizeMin(i) = UseableSize(i)/CHPsize*60;
        OptSize(i) = totalUseSize/OptimalBatSize;
    elseif strcmp(sizeChange,'editSizeMin') 
        UseableSize(i) = SizeMin(i)/60*CHPsize;
        Size(i) = UseableSize(i)/((MaxDOD(i)/100)*(Voltage(i)-DischVoltLoss(i))/Voltage(i));
        OptSize(i) = totalSizeMin/OptimalBatSize;
    elseif strcmp(sizeChange,'OptSize')
        UseableSize(i) = OptSize*OptimalBatSize*(Size(i)/totalSize);
        SizeMin(i) = UseableSize(i)/CHPsize*60;  
        Size(i) = UseableSize(i)/((MaxDOD(i)/100)*(Voltage(i)-DischVoltLoss(i))/Voltage(i));
    end
    Battery(i).OptSize = OptSize(i)*100; 
    Battery(i).Size = Size(i); 
    Battery(i).SizeMin = SizeMin(i); 
end
Project.System.Battery = Battery;