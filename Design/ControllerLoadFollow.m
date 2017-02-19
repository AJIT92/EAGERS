function [DGpower]=ControllerLoadFollow(DemandE,SysSize,RampRate,TurnDown)

steps = length(DemandE);
Ts = 8760/steps;

%% Meet Electric Demand
CHPsize=sum(SysSize);
MinPower = sum(SysSize./TurnDown);
DGpower = zeros(steps+1,1);
DGpower(1) = min(CHPsize,DemandE(1));
RampTs = sum(RampRate.*SysSize);
Hz = min(linspace(8760/Ts-1,0,8760/Ts)',ceil(DemandE/RampTs));
for i = 1:1:steps
    if DemandE(i) > DGpower(i)
        nextPow = min([CHPsize DGpower(i)+RampTs DemandE(i)]);
    else nextPow = max([0 DGpower(i)-RampTs DemandE(i)]);
    end
    [Y, I] = min(DemandE(i:i+Hz(i))-(DGpower(i)-(0:1:Hz(i))'*RampTs));
    if Y<0
        DGpower(i+1) = max(0,DGpower(i)+Y);
        for j = 1:1:I-1
            i = i+1;
            DGpower(i+1) = max(MinPower,DGpower(i)-RampTs);
        end
    else DGpower(i+1) = nextPow;
    end

end
DGpower = DGpower(2:end);