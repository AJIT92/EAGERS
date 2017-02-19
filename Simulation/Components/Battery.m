function Out = Battery(t,Y, Inlet,block,string1)
%1 inlets: V
%6 outlets: IbatMax, IbatMin, VbatMax, VbatMin, E, EMax
%2 states: Temp, SOC
Y = Y.*block.Scale;
Vsoc = block.BatMaxV*(block.k1 + block.k2*log(Y(2)) + block.k3*log(1-Y(2)));%open circuit battery voltage
V = Inlet.V - Vsoc;
I = V/block.Resistance;

QConv = (Y(1) - block.Tamb)*block.AmbConvC*block.SurfA/1000;
QRad = ((Y(1))^4 - (block.Tamb)^4)*block.Epsilon*block.Sigma*block.SurfA/1000;


if strcmp(string1,'dY')
    if I >= 0
        eta = 1;%100% discharge efficiency
        if Y(2) > 0%apply leakage current
            Itotal = I + block.ILeak;
        end
    else
        eta = block.ChargeEff;%configurable charge efficiency
        if Y(2) > 0%apply leakage current
            Itotal = I - block.ILeak;
        end
    end
    
    QElec = Itotal*V/1000;
    if I < 0  && Y(2) > 0%leakage current decreases effectiveness of charging current
        Itotal = Itotal +2*block.ILeak;
    end
    
    dY(1) = (QElec - QConv - QRad)/(block.Mass*block.SpecHeat);
    dY(2) = -eta*Itotal/block.Capacity;
    Out = dY./block.Scale;
    
elseif strcmp(string1,'Outlet')
    %% calculate voltage and current restrictions
    delE = block.Mass*block.SpecHeat*(block.TMax - Y(1))/block.tscale;%Remaining Heat Sink energy in kJ divided over block.tscale seconds;
    IMaxT = (1000*(QConv + QRad + delE)/block.Resistance)^0.5;%Positive Current(amps) that would hit TMax after block.tscale seconds;
    
    IMax = min(block.IMax, IMaxT) - block.ILeak;%apply current and temperature limitations
    
    VMax = min((block.BatMaxV*Y(2)+IMax*block.Resistance), block.VMax);%apply voltage limitations
    VMin = max(block.VMin,(block.BatMaxV*Y(2)-IMax*block.Resistance));
    
    IMaxV = (block.BatMaxV*Y(2) - VMax)/block.Resistance;%current at max V
    IMinV = (block.BatMaxV*Y(2) - VMin)/block.Resistance;%current at min V
    
    %% calculate energy in battery
    
    nSoC = linspace(block.SoCMin, Y(2), block.n);
    E = 0;
    for i = 1:(block.n - 1)%numeric integral of the energy in a battery
        ndt = block.Capacity*(nSoC(i+1) - nSOC(i));%differential time is found by dividing this by the current.
        nV1 = block.BatMaxV*(block.k1 + block.k2*log(nSoC(i)) + block.k3*log(1-nSoC(i)));
        nV2 = block.BatMaxV*(block.k1 + block.k2*log(nSoC(i+1)) + block.k3*log(1-nSoC(i+1)));
        nV = (nV1+nV2)/2;%evaluates each section at midpoint
        E = E + nV*ndt;%E = V*I*dt, but ndt = dt/I
    end
    
    nSoC = linspace(block.SoCMin, block.SoCMax, block.n);
    EMax = 0;
    for i = 1:(block.n - 1)    
        ndt = block.Capacity*(nSoC(i+1) - nSOC(i));
        nV1 = block.BatMaxV*(block.k1 + block.k2*log(nSoC(i)) + block.k3*log(1-nSoC(i)));
        nV2 = block.BatMaxV*(block.k1 + block.k2*log(nSoC(i+1)) + block.k3*log(1-nSoC(i+1)));
        nV = (nV1+nV2)/2;
        EMax = EMax + nV*ndt;
    end
    
    %%
    
    Out.IbatMax = IMaxV;%Current from max circuit voltage(charge)
    Out.IbatMin = IMinV;%Current from min circuit voltage(discharge)
    Out.VbatMax = VMax;%max circuit voltage
    Out.VbatMin = VMin;%min circuit voltage
    Out.E = E;
    Out.EMax = EMax;
    
end

