function Out = HeatExchanger(t,Y, Inlet,block,string1)
%Nodal Heat Exchanger model with 3 states per node: Temperature hot, Plate temperature, temperature cold
% Five (4) inlets: {'Cold Flow','Hot Flow','Cold Pout','Hot Pout'}
global Ru Tags

nodes = block.nodes;
%% scale states to physical values
Y = Y.*block.Scale;
ColdFlow = NetFlow(Inlet.ColdIn);
HotFlow = NetFlow(Inlet.HotIn);

%% seperate out temperatures
Tcold = Y(1:nodes);
Tplate = Y(nodes+1:2*nodes);
Thot = Y(2*nodes+1:3*nodes);
PcoldIn = Y(3*nodes+1);
PhotIn = Y(3*nodes+2);

NcoldOut = block.PfactorCold*(PcoldIn-Inlet.ColdPout);%total cold flow out
NhotOut = block.PfactorHot*(PhotIn-Inlet.HotPout);%total cold flow out

specCold = fieldnames(Inlet.ColdIn);
specHot = fieldnames(Inlet.HotIn);

if strcmp(string1,'Outlet')
    %% Outlet Ports
    Out.ColdOut.T  = mean(Tcold(block.ColdFlowDir(:,end),1));
    for i = 1:1:length(specCold)
        if ~strcmp(specCold{i},'T')
            Out.ColdOut.(specCold{i}) = Inlet.ColdIn.(specCold{i})*NcoldOut/ColdFlow;
        end
    end

    Out.HotOut.T  = mean(Thot(block.HotFlowDir(:,end),1));
    for i = 1:1:length(specHot)
        if ~strcmp(specHot{i},'T')
            Out.HotOut.(specHot{i}) = Inlet.HotIn.(specHot{i})*NhotOut/HotFlow;
        end
    end
    Out.ColdPin = PcoldIn;
    Out.HotPin = PhotIn;
    Tags.(block.name).ColdOut = Tcold(block.ColdFlowDir(:,end),1);
    Tags.(block.name).HotOut = Thot(block.HotFlowDir(:,end),1);
    %% calculate effectiveness & imbalance
    ColdMax = Inlet.ColdIn;
    ColdMax.T = Inlet.HotIn.T;
    HotMin = Inlet.HotIn;
    HotMin.T = Inlet.ColdIn.T;
    HinCold = enthalpy(Inlet.ColdIn);
    HoutCold = enthalpy(Out.ColdOut);
    HinHot = enthalpy(Inlet.HotIn);
    HoutHot = enthalpy(Out.HotOut);
    QT = HoutCold - HinCold;
    maxQT1 = enthalpy(ColdMax) - HinCold;
    maxQT2 = HinHot - enthalpy(HotMin);
    Tags.(block.name).Effectiveness = QT/min(maxQT1,maxQT2);
    Tags.(block.name).NetImbalance = (HinCold + HinHot) - (HoutCold + HoutHot);
elseif strcmp(string1,'dY')
    %% Cold flow
    ColdOutlet.T = Tcold;
    for j = 1:1:length(block.ColdFlowDir(1,:));%1:columns
        k = block.ColdFlowDir(:,j);
        r = length(k);
        if j==1
            ColdInlet.T(k,1) = Inlet.ColdIn.T;
        else
            ColdInlet.T(k,1) = ColdOutlet.T(kprev);
        end
        for i = 1:1:length(specCold)
            if ~strcmp(specCold{i},'T')
                ColdInlet.(specCold{i})(k,1) = Inlet.ColdIn.(specCold{i})/r;
                ColdOutlet.(specCold{i})(k,1) = Inlet.ColdIn.(specCold{i})/r;
            end
        end
        kprev = k;
    end

    %% Hot flow
    HotOutlet.T = Thot;
    for j = 1:1:length(block.HotFlowDir(1,:));%1:columns
        k = block.HotFlowDir(:,j);
        r = length(k);
        if j==1
            HotInlet.T(k,1) = Inlet.HotIn.T;
        else
            HotInlet.T(k,1) = HotOutlet.T(kprev);
        end
        for i = 1:1:length(specHot)
            if ~strcmp(specHot{i},'T')
                HotInlet.(specHot{i})(k,1) = Inlet.HotIn.(specHot{i})/r;
                HotOutlet.(specHot{i})(k,1) = Inlet.HotIn.(specHot{i})/r;
            end
        end
        kprev = k;
    end
    QT = block.HTconv*Y(1:3*nodes) + block.HTcond*Y(1:3*nodes);
    dY = 0*Y;
    %energy flows & sepcific heats
    [~,HoutCold] = enthalpy(ColdOutlet);
    [~,HinCold] = enthalpy(ColdInlet);
    Cp_cold = SpecHeat(ColdOutlet);
    [~,HoutHot] = enthalpy(HotOutlet);
    [~,HinHot] = enthalpy(HotInlet);
    Cp_hot = SpecHeat(HotOutlet);
    
    % time constants for states
    tC1 = (block.Vol_Cold*Cp_cold(k)*PcoldIn./(Ru*Tcold(k)));
    tC2 = (block.Mass*block.Solid_SpecHeat);
    tC3 = (block.Vol_Hot*Cp_hot(k)*PhotIn./(Ru*Thot(k)));
    
    for i=1:1:length(block.ColdFlowDir(1,:))
        k = block.ColdFlowDir(:,i);
        dY(k)= (QT(k) + HinCold(k) - HoutCold(k))./tC1; %Cold flow
        if i>1
            dY(k) = dY(k)+dY(kprev);
        end
        kprev = k;
    end
    dY(nodes+1:2*nodes)= QT(nodes+1:2*nodes)./tC2;  % Solid
    for i=1:1:length(block.HotFlowDir(1,:))
        k = block.HotFlowDir(:,i);
        dY(2*nodes+k)= (QT(2*nodes+k) + HinHot(k) - HoutHot(k))./tC3; %Hot flow
        if i>1
            dY(2*nodes+k) = dY(2*nodes+k)+dY(2*nodes+kprev);
        end
        kprev = k;
    end
    n = 3*nodes;
    %% Pressure
    dY(n+1) = (ColdFlow-NcoldOut)*Ru*Inlet.ColdIn.T/(block.Vol_Cold);%working with total flow rates 
    dY(n+2) = (HotFlow-NhotOut)*Ru*Inlet.HotIn.T/(block.Vol_Hot);%working with total flow rates 

    dY = dY./block.Scale;
    Out = dY;
end