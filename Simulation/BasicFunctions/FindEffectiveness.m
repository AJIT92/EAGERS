function [Effectiveness,Imbalance] = FindEffectiveness(Cold,Hot,ColdOut,HotOut)
%Given the inlet streams (cold & hot, and one of the outlet streams (ColdOut or HotOut), calculate the effectiveness
% if both ColdOut and HotOut are provided, calculate the net energy imbalance

ColdMax = Cold;
ColdMax.T = Hot.T;
HotMin = Hot;
HotMin.T = Cold.T;
maxQT1 = enthalpy(ColdMax) - enthalpy(Cold);
maxQT2 = enthalpy(Hot) - enthalpy(HotMin);
if ~isempty(ColdOut)
    scaleCold = NetFlow(ColdOut)/NetFlow(Cold);
    QT = enthalpy(ColdOut)*scaleCold - enthalpy(Cold);
    if ~isempty(HotOut)
        scaleHot = NetFlow(HotOut)/NetFlow(Hot);
        Imbalance = (enthalpy(Cold) + enthalpy(Hot)) - (enthalpy(ColdOut)*scaleCold + enthalpy(HotOut)*scaleHot);
    else
        Imbalance = 0;
    end
else
    scaleHot = NetFlow(HotOut)/NetFlow(Hot);
    QT = enthalpy(Hot) - enthalpy(HotOut)*scaleHot;
end
Effectiveness = QT/min(maxQT1,maxQT2);