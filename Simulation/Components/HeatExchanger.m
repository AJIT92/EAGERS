function Out = HeatExchanger(t,Y, Inlet,block,string1)
%Nodal Heat Exchanger model with 3 states per node: Temperature hot, Plate temperature, temperature cold
% Five (4) inlets: {'Cold Flow','Hot Flow','Cold Pout','Hot Pout'}
global Ru Tags
ColdFlow = NetFlow(Inlet.Flow1);
HotFlow = NetFlow(Inlet.Flow2);

%% seperate out temperatures
nodes = block.nodes;
Tcold = Y(1:nodes);
Tplate = Y(nodes+1:2*nodes);
Thot = Y(2*nodes+1:3*nodes);
PcoldIn = Y(3*nodes+1);
PhotIn = Y(3*nodes+2);

NcoldOut = block.PfactorCold*(PcoldIn-Inlet.ColdPout);%total cold flow out
NhotOut = block.PfactorHot*(PhotIn-Inlet.HotPout);%total cold flow out

specCold = fieldnames(Inlet.Flow1);
specHot = fieldnames(Inlet.Flow2);

if strcmp(string1,'Outlet')
    ColdOut.T  = mean(Tcold(block.Flow1Dir(:,end),1));
    for i = 1:1:length(specCold)
        if ~strcmp(specCold{i},'T')
            ColdOut.(specCold{i}) = Inlet.Flow1.(specCold{i})*NcoldOut/ColdFlow;
        end
    end
    HotOut.T  = mean(Thot(block.Flow2Dir(:,end),1));
    for i = 1:1:length(specHot)
        if ~strcmp(specHot{i},'T')
            HotOut.(specHot{i}) = Inlet.Flow2.(specHot{i})*NhotOut/HotFlow;
        end
    end
    %% Outlet Ports
    Out.ColdOut = ColdOut;
    Out.HotOut = HotOut;
    
    Out.ColdPin = PcoldIn;
    Out.HotPin = PhotIn;
    Tags.(block.name).ColdOut = Tcold(block.Flow1Dir(:,end),1);
    Tags.(block.name).HotOut = Thot(block.Flow2Dir(:,end),1);
    [Tags.(block.name).Effectiveness , Tags.(block.name).NetImbalance] = FindEffectiveness(Inlet.Flow1,Inlet.Flow2,Out.ColdOut,Out.HotOut); %% calculate effectiveness & imbalance
elseif strcmp(string1,'dY')
    %% Cold flow
    ColdOutlet.T = Tcold;
    for j = 1:1:length(block.Flow1Dir(1,:));%1:columns
        k = block.Flow1Dir(:,j);
        r = length(k);
        if j==1
            ColdInlet.T(k,1) = Inlet.Flow1.T;
        else
            ColdInlet.T(k,1) = ColdOutlet.T(kprev);
        end
        for i = 1:1:length(specCold)
            if ~strcmp(specCold{i},'T')
                ColdInlet.(specCold{i})(k,1) = Inlet.Flow1.(specCold{i})/r;
                ColdOutlet.(specCold{i})(k,1) = Inlet.Flow1.(specCold{i})/r;
            end
        end
        kprev = k;
    end

    %% Hot flow
    HotOutlet.T = Thot;
    for j = 1:1:length(block.Flow2Dir(1,:));%1:columns
        k = block.Flow2Dir(:,j);
        r = length(k);
        if j==1
            HotInlet.T(k,1) = Inlet.Flow2.T;
        else
            HotInlet.T(k,1) = HotOutlet.T(kprev);
        end
        for i = 1:1:length(specHot)
            if ~strcmp(specHot{i},'T')
                HotInlet.(specHot{i})(k,1) = Inlet.Flow2.(specHot{i})/r;
                HotOutlet.(specHot{i})(k,1) = Inlet.Flow2.(specHot{i})/r;
            end
        end
        kprev = k;
    end
    QT = block.HTconv*Y(1:3*nodes) + block.HTcond*Y(1:3*nodes);
    dY = 0*Y;
    %energy flows & sepcific heats
    HoutCold = enthalpy(ColdOutlet);
    HinCold = enthalpy(ColdInlet);
    Cp_cold = SpecHeat(ColdOutlet);
    HoutHot = enthalpy(HotOutlet);
    HinHot = enthalpy(HotInlet);
    Cp_hot = SpecHeat(HotOutlet);
    
    % time constants for states
    tC1 = (block.Vol_Cold*Cp_cold(k)*PcoldIn./(Ru*Tcold(k)));
    tC2 = (block.Mass*block.Solid_SpecHeat);
    tC3 = (block.Vol_Hot*Cp_hot(k)*PhotIn./(Ru*Thot(k)));
    
    for i=1:1:length(block.Flow1Dir(1,:))
        k = block.Flow1Dir(:,i);
        dY(k)= (QT(k) + HinCold(k) - HoutCold(k))./tC1; %Cold flow
        if i>1
            dY(k) = dY(k)+dY(kprev);
        end
        kprev = k;
    end
    dY(nodes+1:2*nodes)= QT(nodes+1:2*nodes)./tC2;  % Solid
    for i=1:1:length(block.Flow2Dir(1,:))
        k = block.Flow2Dir(:,i);
        dY(2*nodes+k)= (QT(2*nodes+k) + HinHot(k) - HoutHot(k))./tC3; %Hot flow
        if i>1
            dY(2*nodes+k) = dY(2*nodes+k)+dY(2*nodes+kprev);
        end
        kprev = k;
    end
    n = 3*nodes;
    %% Pressure
    dY(n+1) = (ColdFlow-NcoldOut)*Ru*Inlet.Flow1.T/(block.Vol_Cold);%working with total flow rates 
    dY(n+2) = (HotFlow-NhotOut)*Ru*Inlet.Flow2.T/(block.Vol_Hot);%working with total flow rates 
    Out = dY;
end