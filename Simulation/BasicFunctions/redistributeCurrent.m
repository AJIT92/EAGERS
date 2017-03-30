function Current = redistributeCurrent(Current,scale,Voltage,localR,SetVoltage)
currentPercOfAvg = Current/sum(Current)*length(Current);
dCurrent = (Voltage-SetVoltage)./localR.*currentPercOfAvg; %change in current to balance voltage
if min(abs(Current./dCurrent))<1
    a = .5*min(abs(Current./dCurrent));%ensure dCurrent is never more than half of a step towards zero
else a = .5;
end
dCurrent = a*dCurrent;
scale2  = (scale*sum(Current))/sum(Current+dCurrent); %change in current to get to new power
Current = (Current+dCurrent)*scale2;%re-distribute current to achieve average voltage, then scale to new total current