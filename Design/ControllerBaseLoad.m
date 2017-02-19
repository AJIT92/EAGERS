function [DGpower]=ControllerBaseLoad(DemandE,SysSize)
steps = length(DemandE);
CHPsize=sum(SysSize);
DGpower(1:steps) = CHPsize;
