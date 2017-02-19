function [DGpower]=ControllerLoadFollow2(DemandE,SysSize,RampRate,TurnDown)

steps = length(DemandE);

%% Meet Electric Demand
CHPsize=sum(SysSize);
MinPower = sum(SysSize./TurnDown);
DGpower = zeros(steps,1);
DGpower(1) = min(CHPsize,DemandE(1));
RampTs = sum(RampRate.*SysSize);
for i = 2:1:steps
    if DemandE(i) > DGpower(i-1)
        DGpower(i) = min([CHPsize DGpower(i-1)+RampTs DemandE(i)]);
    else
        DGpower(i) = max([MinPower DGpower(i-1)-RampTs DemandE(i)]);
    end   
end
