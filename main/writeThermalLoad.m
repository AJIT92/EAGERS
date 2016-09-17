function writeThermalLoad(varargin)
%this function is used at E-hub to simulate a thermal load with the fans
global FanPortWrite DateSim
RealData = GetCurrentData(DateSim);
fwrite(FanPortWrite,num2str(RealData.Demand.H),'char');