function Out = Capacitor(t,Y, Inlet,block,string1)
%1 inlets: V
%6 outlets: IcapMax, IcapMin, VcapMax, VcapMin, E, EMax
%2 states: Temp, Voltage
Y = Y.*block.Scale;
V = Y(2) - Inlet.V;%negative value means charging
I = V/block.Resistance;%negative value means charging


QConv = (Y(1) - block.Tamb)*block.AmbConvC*block.SurfA/1000;
QRad = ((Y(1))^4 - (block.Tamb)^4)*block.Epsilon*block.Sigma*block.SurfA/1000;

if strcmp(string1,'dY')
    dY = 0*Y;
    
    %include effects of current leakage
    Itotal = 0;
    if I >= 0
        if Y(2) > 0
            Itotal = I + block.ILeak;
        end
    else%I < 0
        if Y(2) > 0
            Itotal= I - block.ILeak;
        end
    end
    
    QElec = Itotal*V/1000;%Waste heat in kW
    if I < 0  && Y(2) > 0%leakage current decreases effectiveness of charging current
        Itotal = Itotal +2*block.ILeak;
    end
    
    dY(1) = (QElec - QConv - QRad)/(block.Mass*block.SpecHeat);%temp change
    dY(2) = -Itotal/block.Capacitance;%voltage change
    Out = dY./block.Scale;
elseif strcmp(string1,'Outlet')
    delE = block.Mass*block.SpecHeat*(block.TMax - Y(1))/block.tscale;%Remaining Heat Sink energy in kJ divided over block.tscale seconds;
    IMaxT = (1000*(QConv + QRad + delE)/block.Resistance)^0.5;%Positive Current(amps) that would hit TMax after block.tscale seconds;
    
    IMax = min(block.IMax, IMaxT) - block.ILeak;%apply current and temperature limitations
    
    VMax = min((Y(2)+IMax*block.Resistance), block.VMax);%apply voltage limitations
    VMin = max(0,(Y(2)-IMax*block.Resistance));
    
    IMaxV = (Y(2) - VMax)/block.Resistance;%current at max V
    IMinV = (Y(2) - VMin)/block.Resistance;%current at min V
    
    E = (0.5*block.Capacitance*Y(2)^2)/1000;%Energy stored in capacitor in kJ
    EMax = (0.5*block.Capacitance*block.VMax^2)/1000;
    
    Out.IcapMax = IMaxV;%Current from max circuit voltage(charge)
    Out.IcapMin = IMinV;%Current from min circuit voltage(discharge)
    Out.VcapMax = VMax;%max circuit voltage
    Out.VcapMin = VMin;%min circuit voltage
    Out.E = E;
    Out.EMax = EMax;
end

