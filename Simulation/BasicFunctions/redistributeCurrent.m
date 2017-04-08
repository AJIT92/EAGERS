function [block,error,scale] = redistributeCurrent(block,Inlet,count,firstSolve)
%find the distribution of current that matches the desired power or current, and equalizes voltages
global Tags
Voltage = Tags.(block.name).nVoltage;
OldVoltage = block.Voltage;
OldCurrent = abs(sum(block.Current.H2 + block.Current.CO));
block.Voltage = sum(Tags.(block.name).nVoltage)/block.nodes;
block.Current.CO = Tags.(block.name).I_CO;
block.Current.H2 = Tags.(block.name).I_H2;
netCurrent = block.Current.H2+block.Current.CO;
TotCurrent = abs(sum(netCurrent));
localR = Tags.(block.name).LocalOhmic./abs(netCurrent);

if strcmp(block.Specification,'power density')
    error = (block.RatedStack_kW - Tags.(block.name).Power)/block.RatedStack_kW;
    if count>1
        if firstSolve==1 && error>1e-3
            dP_di = max(.7*(Tags.(block.name).Power - OldCurrent*OldVoltage*block.Cells/1000)/(TotCurrent - OldCurrent),.7*Tags.(block.name).Power/(TotCurrent));%change in power with change in current
        else
            dP_di = max(1.15*((Tags.(block.name).Power - OldCurrent*OldVoltage*block.Cells/1000)/(TotCurrent - OldCurrent)),.7*Tags.(block.name).Power/(TotCurrent));
        end
        scale = 1+ error*block.RatedStack_kW/dP_di/TotCurrent;
    else % first time through
        scale = (block.RatedStack_kW*1000/block.Cells/block.Voltage)/TotCurrent; %total current it should have at this new voltage/ total current specified right now
    end
elseif strcmp(block.Specification,'voltage')
    error = (block.SpecificationValue - block.Voltage)/block.Voltage;
%     if count>1
%         dV_di = -Tags.(block.name).ASR/block.A_Node/100^2;
%         scale = 1 + (block.SpecificationValue - block.Voltage)/dV_di/TotCurrent;
%     else % first time through
%         dV_di = -localR;
%         scale =1+sum((Tags.(block.name).nVoltage - block.SpecificationValue)./localR)/TotCurrent;
%     end
    scale = 1 + 1.5*(block.SpecificationValue - block.Voltage)/mean(-localR)/TotCurrent;
elseif strcmp(block.Specification,'current density')
    scale = abs(Inlet.NetCurrent)/TotCurrent;
    error = (TotCurrent - abs(Inlet.NetCurrent))/TotCurrent;
end
ratio = block.Current.CO./netCurrent;
    
%% start by maintaining same current in each row, then allow row voltages to balance (prevents fuel starvation in a row during initialization)
currentPercOfAvg = netCurrent/(sum(netCurrent)/length(netCurrent));
dCurrent = (Voltage-block.Voltage)./localR.*currentPercOfAvg; %change in current to balance voltage
if min(abs(netCurrent./dCurrent))<1
    a = .5*min(abs(netCurrent./dCurrent));%ensure dCurrent is never more than half of a step towards zero
else a = .5;
end
dCurrent = a*dCurrent;
scale2  = (scale*sum(netCurrent))/sum(netCurrent+dCurrent); %change in current to get to new power
netCurrent = (netCurrent+dCurrent)*scale2;%re-distribute current to achieve average voltage, then scale to new total current

block.Current.CO = ratio.*netCurrent;
block.Current.H2 = (1-ratio).*netCurrent;