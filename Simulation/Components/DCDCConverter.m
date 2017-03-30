function Out = DCDCConverter(t,Y, Inlet,block,string1)
%A very simple DC/DC Converter
%0 States
%6 Inlets: ILoad, VLoad, VMax, Ivmax, VMin, Ivmin
%3 Outlet: Imin, Imax, VSource
if strcmp(string1,'Outlet')
    Ivmax = block.Effic*Inlet.Ivmax*Inlet.VMax/Inlet.VLoad;%DC bus current when v = vmax
    Ivmin = block.Effic*Inlet.Ivmin*Inlet.VMin/Inlet.VLoad;%DC bus current when v = vmin
    
    %if current is negative, efficiency losses are applied at dc bus
    if Inlet.Ivmax < 0
        Ivmax = Inlet.Ivmax*Inlet.VMax/(Inlet.VLoad*block.Effic);
    end
    if Inlet.Ivmin < 0
        Ivmin = Inlet.Ivmin*Inlet.VMin/(Inlet.VLoad*block.Effic);
    end
    
    V = Inlet.VMin +(Inlet.VMax - Inlet.VMin)*(Inlet.ILoad - Ivmin)/(Ivmax - Ivmin);%voltage on source circuit that creates ILoad on DC bus
    
    IMin = min(Ivmax,Ivmin);
    IMax = max(Ivmax,Ivmin);
    
    Out.IMin = IMin;%min current applied to DC bus by source
    Out.IMax = IMax;%max current applied to DC bus by source
    
    Out.VSource = V;%voltage on source circuit
elseif strcmp(string1,'dY')
    %no states
end