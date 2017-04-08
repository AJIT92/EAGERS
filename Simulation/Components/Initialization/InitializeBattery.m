function block = InitializeBattery(varargin)
%A simple capacitor model
%1 inlets: V
%6 outlets: IbatMax, IbatMin, VbatMax, VbatMin, E, EMax
%2 states: Temp, SoC(state of charge)
block = varargin{1};
if length(varargin) ==1 %first initialization
    block.Tamb = 305;%ambient air temp
    block.AmbConvC = 5;%ambient convection coefficient
    block.Epsilon = .8;%radiation epsilon value
    block.Sigma = 5.67e-8;%Radiation sigma value
    block.SpecHeat = .5;%Specific Heat of Capacitor kJ/(kg*K)

    %% can be pulled out of initialization function
    
    block.Mass = 1;
    block.Capacity = 1;%in amp hours
    block.BatMaxV = 1;%in volts
    block.n = 50;% number of integration points in battery energy storage calculation
    block.SoCMin = 0.3;%restrict SoC to maximize battery life
    block.SoCMax = 0.7;%restrict SoC to maximize battery life
    block.Resistance = 1;%in ohms
    block.ChargeEff = 1;%what percentage of charging current charges battery SoC
    block.ILeak = 1;%leakage current in amps
    block.TMax = 338;%in K
    block.IMax = 1;% in amps
    block.InitialSoC = 0.5;
    block.tscale = 1;
    
    %Vsoc = block.BatMaxV*(block.k1 + block.k2*log(Y(2)) + block.k3*log(1-Y(2)));%open circuit battery voltage
    block.k1 = 0.9498;%k constants found by fitting battery data
    block.k2 = 0.03168;
    block.k3 = -0.004618;
    
    %block.VMax = block.BatMaxV*(block.k1 + block.k2*log(block.SoCMax) + block.k3*log(1-block.SoCMax));
    %block.VMin = block.BatMaxV*(block.k1 + block.k2*log(block.SoCMin) + block.k3*log(1-block.SoCMin));
    
    block.VMax = block.BatMaxV;
    block.VMin = 0;
    block.SurfA = 2*block.Width*block.Length + 2* block.Width*block.Height + 2*block.Length*block.Height;%Assume rectangular prism
    
    block.Scale = [block.Tamb, 1];%SoC is between 0 and 1
    block.IC = [1,block.InitialSoC];
    
    block.InletPorts = {'V'};
    block.V.IC = block.BatMaxV*(block.k1 + block.k2*log(block.Scale(2)*block.IC(2)) + block.k3*log(1-block.Scale(2)*block.IC(2)));
    
    block.OutletPorts = {'IbatMax','IbatMin','VbatMax','VbatMin','E','EMax'};
    block.IbatMax.IC = -block.IMax;%current at max applied voltage
    block.IbatMin.IC = block.IMax;%current at min applied voltage
    block.VbatMax.IC = block.Scale(2)*block.IC(2) + (block.IbatMax.IC*block.Resistance);
    block.VbatMin.IC = block.Scale(2)*block.IC(2) + (block.IbatMin.IC*block.Resistance);
    
    nSoC = linspace(block.SoCMin, block.Scale(2)*block.IC(2), block.n);
    E = 0;
    for i = 1:(block.n - 1)%numeric integral of the energy in a battery
        ndt = block.Capacity*(nSoC(i+1) - nSOC(i));%differential time is found by dividing this by the current.
        nV1 = block.BatMaxV*(block.k1 + block.k2*log(nSoC(i)) + block.k3*log(1-nSoC(i)));
        nV2 = block.BatMaxV*(block.k1 + block.k2*log(nSoC(i+1)) + block.k3*log(1-nSoC(i+1)));
        nV = (nV1+nV2)/2;%evaluates each section at midpoint
        E = E + nV*ndt;%E = V*I*dt, but ndt = dt/I
    end
    block.E.IC = E;
     
    nSoC = linspace(block.SoCMin, block.SoCMax, block.n);
    EMax = 0;
    for i = 1:(block.n - 1)    
        ndt = block.Capacity*(nSoC(i+1) - nSOC(i));
        nV1 = block.BatMaxV*(block.k1 + block.k2*log(nSoC(i)) + block.k3*log(1-nSoC(i)));
        nV2 = block.BatMaxV*(block.k1 + block.k2*log(nSoC(i+1)) + block.k3*log(1-nSoC(i+1)));
        nV = (nV1+nV2)/2;
        EMax = EMax + nV*ndt;
    end
    block.EMax.IC = EMax;
elseif length(varargin) ==2 %have inlets
    %Inlet = varargin{2};%states don't change during initialization
    
    Y = block.IC.*block.Scale;
    %Vsoc = block.BatMaxV*(block.k1 + block.k2*log(Y(2)) + block.k3*log(1-Y(2)));%open circuit battery voltage
    %V = Inlet.V - Vsoc;
    %I = V/block.Resistance;

    QConv = (Y(1) - block.Tamb)*block.AmbConvC*block.SurfA/1000;
    QRad = ((Y(1))^4 - (block.Tamb)^4)*block.Epsilon*block.Sigma*block.SurfA/1000;
    
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
    
    %update port initial conditions
    block.IbatMax.IC = IMaxV;%Current from max circuit voltage(charge)
    block.IbatMin.IC = IMinV;%Current from min circuit voltage(discharge)
    block.VbatMax.IC = VMax;%max circuit voltage
    block.VbatMin.IC = VMin;%min circuit voltage
    block.E.IC = E;
    block.EMax.IC = EMax;
    
    % need to update block.Scale
end

