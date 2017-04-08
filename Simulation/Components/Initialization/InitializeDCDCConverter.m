function block = InitializeDCDCConverter( varargin )
%A very simple DC/DC converter model
%0 States
%6 Inlets: ILoad, VLoad, VMax, Ivmax, VMin, Ivmin
%3 Outlet: Imin, Imax, VSource
block = varargin{1};
if length(varargin) ==1
    
    block.Efficiency = 0.90;%can be moved to model file
    
    block.InletPorts = {'Ivmax','VMax','Ivmin','VMin','ILoad','VLoad'};
    block.Ivmax.IC = -1;
    block.VMax.IC = 50;
    block.Ivmin.IC = 1;
    block.VMin.IC = 30;
    block.ILoad.IC = 0;
    block.VLoad.IC = 80;
    
    block.OutletPorts = {'Imin','Imax','VSource'};
    block.Imin.IC = -1;
    block.Imax.IC = 1;
    block.VSource.IC = 40;
elseif length(varargin) ==2
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
    
    block.Imin.IC = IMin;
    block.Imax.IC = IMax;
    block.VSource.IC = V;
    
    % need to update block.Scale
end

