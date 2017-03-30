function Out = Electrolyzer(t,Y, Inlet,block,string1)
%This function models an electorlyzer it takes input Y states and returns dY: 
%Electrolyzer model with many states: Temperatures (oxidizer plate, cathode, electrolyte, anode, fuel plate [reformer]), Cathode species ( [ CO2, H2O], N2, O2) anode species (CH4, CO, CO2, H2, H2O, N2, O2), [Reformer species] [Rate of internal reforming reactions] Current, cathode pressure, anode pressure
% Five (5) inlets: {'NetCurrent','Flow1','Flow2','Flow1Pout','Flow2Pout'}
% Seven (7) outlets: {'Flow1Out','Flow2Out','Flow1Pin','Flow2Pin','MeasureVoltage','MeasureTpen','MeasureTflow1','MeasureTflow2'}
% for an elextrolyzer, Flow 1 is the cathode (steam), and Flow2 is the anode (cooling/heating air)
% Current is negative
global F Ru Tags
nodes = block.nodes;
%% add species that may not be in inlet
inFields = fieldnames(Inlet.Flow1);
for i = 1:1:length(block.F1Spec)
    if ~ismember(block.F1Spec{i},inFields)
        Inlet.Flow1.(block.F1Spec{i}) = 0;
    end
end
inFields = fieldnames(Inlet.Flow2);
for i = 1:1:length(block.F2Spec)
    if ~ismember(block.F2Spec{i},inFields)
        Inlet.Flow2.(block.F2Spec{i}) = 0;
    end
end
%% separate out temperatures
T_plate1 = Y(1:nodes);
Flow1.Outlet.T = Y(nodes+1:2*nodes);
T_Elec = Y(2*nodes+1:3*nodes);
Flow2.Outlet.T = Y(3*nodes+1:4*nodes); 
T_plate2 = Y(4*nodes+1:5*nodes); 
switch block.Reformer
    case 'methanator'
        nT = 6*nodes; % # of temperature states
        Flow3.Outlet.T = Y(5*nodes+1:6*nodes); %only gets used if internal Methanator exists, otherwise these values are actually QT1
    case 'none'
        nT = 5*nodes; % # of temperature states
end
n = nT;

%Current
nCurrent = Y(end-nodes-1:end-2);
%Pressure
P_flow1 = Y(end-1); %pressure
P_flow2 = Y(end); %pressure

%% Cathode
for i = 1:1:length(block.F1Spec)
    Flow1.Outlet.(block.F1Spec{i}) = Y(n+1:n+nodes); n = n+nodes;
end
for j = 1:1:length(block.Flow1Dir(1,:));
    if j==1%first column recieves fresh inlet
        k = block.Flow1Dir(:,1);
        Flow1.Inlet.T(k,1) = Inlet.Flow1.T; 
        for i = 1:1:length(block.F1Spec)
            Flow1.Inlet.(block.F1Spec{i})(k,1) = Inlet.Flow1.(block.F1Spec{i})/block.Cells/length(k); 
        end
    else%subsequent columns recieve outlet of previous column
        k2 = block.Flow1Dir(:,j);
        Flow1.Inlet.T(k2,1) = Flow1.Outlet.T(k);
        for i = 1:1:length(block.F1Spec)
            Flow1.Inlet.(block.F1Spec{i})(k2,1) = Flow1.Outlet.(block.F1Spec{i})(k); 
        end
        k = k2;
    end
end
Flow1Out.T  = mean(Flow1.Outlet.T(block.Flow1Dir(:,end))); %temperature 
for i = 1:1:length(block.F1Spec)
    Flow1Out.(block.F1Spec{i}) = max(0,sum(Flow1.Outlet.(block.F1Spec{i})(block.Flow1Dir(:,end)))*block.Cells);%avoid sending negative outlets
end

%% Anode 
for i = 1:1:length(block.F2Spec)
    Flow2.Outlet.(block.F2Spec{i}) = Y(n+1:n+nodes); n = n+nodes;
end
for j = 1:1:length(block.Flow2Dir(1,:))
    k2 =block.Flow2Dir(:,j);
    if j==1 % first column of fuel flow direction
        Flow2.Inlet.T(k2,1) = Inlet.Flow2.T;
        for i = 1:1:length(block.F2Spec)
            Flow2.Inlet.(block.F2Spec{i})(k2,1) = Inlet.Flow2.(block.F2Spec{i})/block.Cells/length(k2);
        end
    else
        Flow2.Inlet.T(k2,1) = Flow2.Outlet.T(k);
        for i = 1:1:length(block.F2Spec)
            Flow2.Inlet.(block.F2Spec{i})(k2,1) = Flow2.Outlet.(block.F2Spec{i})(k);
        end
    end
    k = k2;
end
Flow2Out.T  = mean(Flow2.Outlet.T(block.Flow2Dir(:,end))); %temperature 
for i = 1:1:length(block.F2Spec)
    Flow2Out.(block.F2Spec{i}) = max(0,sum(Flow2.Outlet.(block.F2Spec{i})(block.Flow2Dir(:,end)))*block.Cells);%avoid sending negative outlets
end

%% Methanator
switch block.Reformer
    case 'methanator' %secondary inlet introduces CO2 stream
        for i = 1:1:length(block.F2Spec)
            Flow3.Outlet.(block.F2Spec{i}) = Y(n+1:n+nodes); n = n+nodes;
        end
        for j = 1:1:length(block.Flow3Dir(1,:))
            if j==1
                k = block.Flow3Dir(:,1);
                Flow3.Inlet.T(k,1) = Inlet.CarbonDioxide.T;
                for i = 1:1:length(block.F2Spec)
                    Flow3.Inlet.(block.F2Spec{i})(k,1) = Flow2.Outlet.(block.F2Spec{i})(block.Flow2Dir(:,1))*block.MethSpacing;
                    if isfield(Inlet.CarbonDioxide,block.F2Spec{i})
                        Flow3.Inlet.(block.F2Spec{i})(k,1) = Flow3.Inlet.(block.F2Spec{i})(k,1) +Inlet.CarbonDioxide.(block.F2Spec{i})/block.Cells/length(k)*block.MethSpacing;
                    end
                end
            else
                k2 = block.Flow3Dir(:,j);
                Flow3.Inlet.T(k2,1) = Flow3.Outlet.T(k);
                for i = 1:1:length(block.F2Spec)
                    Flow3.Inlet.(block.F2Spec{i})(k2,1) = Flow3.Outlet.(block.F2Spec{i})(k);
                end
                k = k2;
            end
        end
end
if strcmp(string1,'Outlet')
    %%Nernst & Losses
    n_an_in = NetFlow(Flow2.Inlet);
    n_an_out = NetFlow(Flow2.Outlet);
    n_cath_in = NetFlow(Flow1.Inlet);
    n_cath_out = NetFlow(Flow1.Outlet);
    AvgX.O2 = (Flow2.Outlet.O2+Flow2.Inlet.O2)./(n_an_in+n_an_out);
    AvgX.H2 = (Flow1.Outlet.H2+Flow1.Inlet.H2)./(n_cath_in+n_cath_out);
    AvgX.H2O = (Flow1.Outlet.H2O+Flow1.Inlet.H2O)./(n_cath_in+n_cath_out);

    FuelCellNernst(nCurrent,T_Elec,P_flow2,AvgX,block);
    Voltage =  sum(Tags.(block.name).nVoltage.*(nCurrent/sum(nCurrent)));
    Tags.(block.name).Voltage = Voltage;
    %% Tags
    H2O_in = Inlet.Flow1.H2O;
    H2O_out = Flow1Out.H2O;
    Tags.(block.name).H2Outilization = (H2O_in - H2O_out)./H2O_in;
    Tags.(block.name).Tpen = T_Elec;
    Tags.(block.name).TcathOut = Flow1Out.T;
    Tags.(block.name).TanodeOut = Flow2Out.T;
    Tags.(block.name).Current = sum(abs(nCurrent));
    Tags.(block.name).StackdeltaT = Flow2Out.T-Inlet.Flow2.T;
    Tags.(block.name).PENavgT = sum(T_Elec)/block.nodes;
    Tags.(block.name).MaxPEN = max(T_Elec);
    Tags.(block.name).PENdeltaT = Tags.(block.name).MaxPEN-min(T_Elec);
    Tags.(block.name).dTdX = (T_Elec-T_Elec(block.HTadjacent(:,2)))/(block.L_Cell/block.columns);
    Tags.(block.name).dTdY = (T_Elec-T_Elec(block.HTadjacent(:,4)))/(block.W_Cell/block.rows);
    Tags.(block.name).MaxdTdX = max(abs([Tags.(block.name).dTdX;Tags.(block.name).dTdY;]));
    %%Outlet Ports
    Out.Flow1Out  = Flow1Out;
    Out.Flow2Out = Flow2Out;
    Out.Flow1Pin = P_flow1;
    Out.Flow2Pin = P_flow2;
    Out.MeasureVoltage = Voltage;
    Out.MeasurePower = Tags.(block.name).Power;
    Out.MeasureTpen = Y(2*nodes+1:3*nodes);
    Out.MeasureTflow1 = Y(nodes+block.Flow1Dir(:,end));
    Out.MeasureTflow2 = Y(3*nodes+block.Flow2Dir(:,end));
elseif strcmp(string1,'dY')  
    Voltage =  Tags.(block.name).Voltage;
    switch block.Reformer
        case 'methanator'
            [~,RefOut] = KineticReformation(block.method,Flow3,P_flow2,block.KineticCoeff3,zeros(block.nodes,1),block);%% Kinetic reaction rates  (WGS is always near equilibrium)
%             R.CH4ref = Rref.CH4;
%             R.WGSref = Rref.WGS;
    end
    
    [h,~] = enthalpy(T_Elec,{'H2','H2O','O2','CO2'});
    Power = Voltage*nCurrent/1000; %cell power in kW
    Qgen = nCurrent/(2000*F).*(h.H2+.5*h.O2-h.H2O)-Power;%kW of heat generated by electrochemistry (per node & per cell)
    switch block.FCtype%ion transport across membrane (total enthalpy)
        case {'SOFC';'SOEC'}
            Qion = nCurrent/(4000*F).*h.O2; %O2 ion crossing over (kW)
        case {'MCFC';'MCEC'}
            Qion = nCurrent/(4000*F).*h.O2 + nCurrent/(2000*F).*h.CO2;% O2 & CO2 ion crossing over
    end
    
    %% Q %% Heat transfer & Generation
    switch block.Reformer
        case 'methanator'
            QT = block.HTcond*Y(1:6*nodes) + block.HTconv*Y(1:6*nodes);
        case {'none'}
            QT = block.HTcond*Y(1:5*nodes) + block.HTconv*Y(1:5*nodes);
    end
    Tags.(block.name).Q_gen = sum(Qgen*block.Cells); %kW of heat generated by electrochemistry
    %energy flows & sepcific heats
    HoutCath = enthalpy(Flow1.Outlet);
    HinCath = enthalpy(Flow1.Inlet);
    HoutAnode = enthalpy(Flow2.Outlet);
    HinAnode = enthalpy(Flow2.Inlet);

    %% %% solve for dY in order of states
    dY = 0*Y;
    %%Temperatures
    dY(1:nodes)= QT(1:nodes)./block.tC(1:nodes);  %Cathode Plate
    for i=1:1:length(block.Flow1Dir(1,:)) %having the downstream nodes change temperature with the upstream nodes prevents propogation issues when taking larger time steps
        k = block.Flow1Dir(:,i);
        dY(nodes+k)= (QT(nodes+k) + HinCath(k) - HoutCath(k)  - Power(k) - Qgen(k) + Qion(k))./block.tC(nodes+k); %Cathode
        if i>1
            dY(nodes+k) = dY(nodes+k)+dY(nodes+kprev);
        end
        kprev = k;
    end
    dY(1+2*nodes:3*nodes)= (QT(1+2*nodes:3*nodes) + Qgen)./block.tC(2*nodes+1:3*nodes); %Electrolyte Plate
    for i=1:1:length(block.Flow2Dir(1,:))
        k = block.Flow2Dir(:,i);
        dY(3*nodes+k)= (QT(3*nodes+k) + HinAnode(k) - HoutAnode(k) - Qion(k))./block.tC(3*nodes+k);  %Anode
        if i>1
            dY(3*nodes+k) = dY(3*nodes+k)+dY(3*nodes+kprev);
        end
        kprev = k;
    end
    dY(1+4*nodes:5*nodes)= QT(1+4*nodes:5*nodes)./block.tC(4*nodes+1:5*nodes);  %Anode Plate

    n =nT;
    %%Cathode Species
    for i = 1:1:length(block.F1Spec)
        if strcmp(block.F1Spec{i},'H2O')
            dY(n+1:n+nodes)= (Flow1.Inlet.H2O - Flow1.Outlet.H2O + nCurrent/(2*F*1000))./block.tC(n+1:n+nodes);  %H2O species concentration with CO2 crossover
        elseif strcmp(block.F1Spec{i},'H2') 
            dY(n+1:n+nodes)= (Flow1.Inlet.H2 - Flow1.Outlet.H2 - nCurrent/(2*F*1000))./block.tC(n+1:n+nodes);%H2 species concentration with O2 crossover
        else
            dY(n+1:n+nodes)= (Flow1.Inlet.(block.F1Spec{i}) - Flow1.Outlet.(block.F1Spec{i}))./block.tC(n+1:n+nodes);%all other species concentration
        end
        n = n+nodes;
    end

    %%Anode Species
    for i = 1:1:length(block.F2Spec)
        if strcmp(block.F2Spec{i},'O2')
            dY(n+1:n+nodes)= (Flow2.Inlet.O2 - Flow2.Outlet.O2 - nCurrent/(4*F*1000))./block.tC(n+1:n+nodes);%O2 species concentration with O2 crossover
        else
            dY(n+1:n+nodes)= (Flow2.Inlet.(block.F2Spec{i}) - Flow2.Outlet.(block.F2Spec{i}))./block.tC(n+1:n+nodes); %all species concentration
        end
        n = n+nodes; 
    end
   
    %%Methanator
    switch block.Reformer
        case 'methanator'
            [~,HoutReform] = enthalpy(Flow3.Outlet);
            [~,HinReform] = enthalpy(Flow3.Inlet);
            for i=1:1:length(block.ReformFlowDir(1,:))
                k = block.ReformFlowDir(:,i);
                dY(5*nodes+k)= (block.RefSpacing*QT(5*nodes+k) + HinReform(k) - HoutReform(k))./block.tC(5*nodes+k);  %Fuel Methanator Channels
                if i>1
                    dY(5*nodes+k) = dY(5*nodes+k)+dY(5*nodes+kprev);
                end
                kprev = k;
            end
            for i = 1:1:length(block.F2Spec)
                dY(n+1:n+nodes)= (RefOut.(block.F2Spec{i}) - Flow3.Outlet.(block.F2Spec{i}))./block.tC(n+1:n+nodes);   %all species concentrations
                n = n+nodes;
            end
    end
    %%Current % note sign convention of current is negative for electrolyzers !!
    dY(n+1:n+nodes) = (Inlet.NetCurrent - sum(nCurrent) + (Tags.(block.name).nVoltage - Voltage)./Tags.(block.name).ASR.*(block.A_Cell*100^2))./block.tC(end-nodes-1:end-2); n = n+nodes; %error in A/cm^2 * area  
    Power = Voltage*abs(Inlet.NetCurrent)*block.Cells/1000;
    Efficiency = 240424*Flow1Out.H2/Power;
    Tags.(block.name).Efficiency = Efficiency;

    %%Pressure
    Nanode = block.PfactorAnode*max(0.01,(P_flow2-Inlet.Flow2Pout));%total anode flow out
    Ncath = block.PfactorCath*max(0.1,(P_flow1-Inlet.Flow1Pout));%total cathode flow out
    dY(n+1) = (NetFlow(Inlet.Flow1)-Ncath)*Ru*Inlet.Flow1.T/block.tC(n+1);%working with total flow rates so must multiply by nodes & cells
    dY(n+2) = (NetFlow(Inlet.Flow2)-Nanode)*Ru*Inlet.Flow2.T/block.tC(n+2);
    Out = dY;
end