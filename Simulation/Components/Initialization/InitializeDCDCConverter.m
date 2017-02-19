function block = InitializeDCDCConverter( varargin )
%A very simple DC/DC converter model
%0 States
%6 Inlets: ILoad, VLoad, VMax, Ivmax, VMin, Ivmin
%3 Outlet: Imin, Imax, VSource
block = varargin{1};
if length(varargin) ==1
    
    block.Efficiency = 0.90;%can be moved to model file
    
    block.PortNames = {'Ivmax','VMax','Ivmin','VMin','ILoad','VLoad','Imin','Imax','VSource'};
    
    block.Ivmax.type = 'in';
    block.Ivmax.IC = -1;
    
    block.VMax.type = 'in';
    block.VMax.IC = 50;
    
    block.Ivmin.type = 'in';
    block.Ivmin.IC = 1;
    
    block.VMin.type = 'in';
    block.VMin.IC = 30;
    
    block.ILoad.type = 'in';
    block.ILoad.IC = 0;
    
    block.VLoad.type = 'in';
    block.VLoad.IC = 80;
    
    block.Imin.type = 'out';
    block.Imin.IC = -1;
    
    block.Imax.type = 'out';
    block.Imax.IC = 1;
    
    block.VSource.type = 'out';
    block.VSource.IC = 40;
    
    for i = 1:1:length(block.PortNames)
        if length(block.connections)<i || isempty(block.connections{i})
            block.(block.PortNames{i}).connected={};
        else
            if ischar(block.connections{i})
                block.(block.PortNames{i}).connected = block.connections(i);
            else
                block.(block.PortNames{i}).IC = block.connections{i};
                block.(block.PortNames{i}).connected={};
            end
        end
    end
    
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

