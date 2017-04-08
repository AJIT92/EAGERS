function Out = FuelCell(t,Y, Inlet,block,string1)
%This function models a fuel cell  it takes input Y states and returns dY: 
%FC model with many states: Temperatures (oxidizer plate, cathode, electrolyte, anode, fuel plate [reformer]), Cathode species ( [ CO2, H2O], N2, O2) anode species (CH4, CO, CO2, H2, H2O, N2, O2), [Reformer species] [Rate of internal reforming reactions] Current, cathode pressure, anode pressure
% Five (5) inlets: {'NetCurrent','Flow1','Flow2','Flow2Pout','Flow1Pout'}
% Seven (7) outlets: {'Flow2Out','Flow1Out','Flow2Pin','Flow1Pin','MeasureCurrent','MeasureTpen','MeasureTflow1','MeasureTflow2'}
% Fuel cell current is positive, Flow 1 is the anode (fuel), and Flow2 is the cathode cooling, flow 3 is the internal reformer (upstream of the anode)
% Electrolyzer current is negative, Flow 1 is the steam input side, and Flow2 is the cooling/heating air, flow 3 is the internal methanation (CO2 injection downstream of the anode)
global F Ru Tags
%% add species that may not be in inlet
inFields = fieldnames(Inlet.Flow1);
for i = 1:1:length(block.Spec1)
    if ~ismember(block.Spec1{i},inFields)
        Inlet.Flow1.(block.Spec1{i}) = 0;
    end
end
inFields = fieldnames(Inlet.Flow2);
for i = 1:1:length(block.Spec2)
    if ~ismember(block.Spec2{i},inFields)
        Inlet.Flow2.(block.Spec2{i}) = 0;
    end
end

nodes = block.nodes;
%% seperate out temperatures
T_plate1 = Y(1:nodes);
Flow2.Outlet.T = Y(nodes+1:2*nodes);
T_Elec = Y(2*nodes+1:3*nodes);
Flow1.Outlet.T = Y(3*nodes+1:4*nodes); 
T_plate2 = Y(4*nodes+1:5*nodes); 
switch block.Reformer
    case 'internal'
        nT = 6*nodes; % # of temperature states
        Flow3.Outlet.T = Y(5*nodes+1:6*nodes); %only gets used if internal reformer exists, otherwise these values are actually QT1
    case 'direct';
        nT = 5*nodes; % # of temperature states
end
n = nT;

%Current
nCurrent = Y(end-nodes-1:end-2);
%Pressure
P_flow1 = Y(end-1); %pressure
P_flow2 = Y(end); %pressure

%% Cathode
for i = 1:1:length(block.Spec2)
    Flow2.Outlet.(block.Spec2{i}) = Y(n+1:n+nodes); n = n+nodes;
end
for j = 1:1:length(block.Flow2Dir(1,:));
    if j==1%first column recieves fresh inlet
        k = block.Flow2Dir(:,1);
        Flow2.Inlet.T(k,1) = Inlet.Flow2.T; 
        for i = 1:1:length(block.Spec2)
            Flow2.Inlet.(block.Spec2{i})(k,1) = Inlet.Flow2.(block.Spec2{i})/block.Cells/length(k); 
        end
    else%subsequent columns recieve outlet of previous column
        k2 = block.Flow2Dir(:,j);
        Flow2.Inlet.T(k2,1) = Flow2.Outlet.T(k);
        for i = 1:1:length(block.Spec2)
            Flow2.Inlet.(block.Spec2{i})(k2,1) = Flow2.Outlet.(block.Spec2{i})(k); 
        end
        k = k2;
    end
    if block.ClosedCathode %closed end cathode
        Flow2.Outlet.O2(k,1) = Flow2.Inlet.O2(k,1) - nCurrent(k)/(4*F*1000); 
        if any(strcmp(block.FCtype,{'MCFC','MCEC'}))
            Flow2.Outlet.CO2(k,1) = Flow2.Inlet.CO2(k,1) - nCurrent(k)/(2*F*1000); 
        end
    end
end

%% Anode 
for i = 1:1:length(block.Spec1)
    Flow1.Outlet.(block.Spec1{i}) = Y(n+1:n+nodes); n = n+nodes;
end
%% Internal reformer
switch block.Reformer
    case 'internal'
        for i = 1:1:length(block.Spec1)
            Flow3.Outlet.(block.Spec1{i}) = Y(n+1:n+nodes); n = n+nodes;
        end
end

Flow1Out.T  = mean(Flow1.Outlet.T(block.Flow1Dir(:,end))); %temperature 
for i = 1:1:length(block.Spec1)
    Flow1Out.(block.Spec1{i}) = max(0,sum(Flow1.Outlet.(block.Spec1{i})(block.Flow1Dir(:,end)))*block.Cells);%avoid sending negative outlets
end
Flow2Out.T  = mean(Flow2.Outlet.T(block.Flow2Dir(:,end))); %temperature 
for i = 1:1:length(block.Spec2)
    Flow2Out.(block.Spec2{i}) = max(0,sum(Flow2.Outlet.(block.Spec2{i})(block.Flow2Dir(:,end)))*block.Cells);%avoid sending negative outlets
end
%% Reformer
switch block.Reformer
    case 'internal'
        for j = 1:1:length(block.Flow3Dir(1,:))
            if j==1
                k = block.Flow3Dir(:,1);
                Flow3.Inlet.T(k,1) = Inlet.Flow1.T;
                for i = 1:1:length(block.Spec1)
                    Flow3.Inlet.(block.Spec1{i})(k,1) = Inlet.Flow1.(block.Spec1{i})/block.Cells/length(k)*block.RefSpacing;
                end
            else
                k2 = block.Flow3Dir(:,j);
                Flow3.Inlet.T(k2,1) = Flow3.Outlet.T(k);
                for i = 1:1:length(block.Spec1)
                    Flow3.Inlet.(block.Spec1{i})(k2,1) = Flow3.Outlet.(block.Spec1{i})(k);
                end
                k = k2;
            end
        end
end
%% Anode
for j = 1:1:length(block.Flow1Dir(1,:))
    k2 =block.Flow1Dir(:,j);
    if j==1 % first column of fuel flow direction
        switch block.Reformer
            case 'internal'
                Flow1.Inlet.T(k2,1) =  Flow3.Outlet.T(k);
            case {'direct';'adiabatic'}
                Flow1.Inlet.T(k2,1) = Inlet.Flow1.T;
        end
        for i = 1:1:length(block.Spec1)
            switch block.Reformer
                case 'internal'
                    Flow1.Inlet.(block.Spec1{i})(k2,1) = Flow3.Outlet.(block.Spec1{i})(k)/block.RefSpacing;%Species flows coming into anode
                case {'direct';'adiabatic'}
                    Flow1.Inlet.(block.Spec1{i})(k2,1) = Inlet.Flow1.(block.Spec1{i})/block.Cells/length(k2);
            end
        end
    else
        Flow1.Inlet.T(k2,1) = Flow1.Outlet.T(k);
        for i = 1:1:length(block.Spec1)
            Flow1.Inlet.(block.Spec1{i})(k2,1) = Flow1.Outlet.(block.Spec1{i})(k);
        end
    end
    k = k2;
end

if strcmp(string1,'Outlet')
    %% Nernst & Losses
    FuelCellNernst(Flow1,Flow2,nCurrent,T_Elec,P_flow2,block)
    Tags.(block.name).Voltage = sum(Tags.(block.name).nVoltage.*(nCurrent/sum(nCurrent)));
    %% Tags
    H2_in = sum(4*Flow1.Inlet.CH4(block.Flow1Dir(:,1))+ Flow1.Inlet.CO(block.Flow1Dir(:,1)) + Flow1.Inlet.H2(block.Flow1Dir(:,1)));
    H2_out = sum(4*Flow1.Outlet.CH4(block.Flow1Dir(:,end))+ Flow1.Outlet.CO(block.Flow1Dir(:,end)) + Flow1.Outlet.H2(block.Flow1Dir(:,end)));
    Tags.(block.name).H2utilization = (H2_in - H2_out)./H2_in;
    Tags.(block.name).O2utilization = sum(Flow2.Inlet.O2(block.Flow2Dir(:,1)) - Flow2.Outlet.O2(block.Flow2Dir(:,end)))/sum(Flow2.Inlet.O2(block.Flow2Dir(:,1)));
    if strcmp(block.FCtype,'MCFC')
        Tags.(block.name).CO2utilization = sum(Flow2.Inlet.CO2(block.Flow2Dir(:,1)) - Flow2.Outlet.CO2(block.Flow2Dir(:,end)))/sum(Flow2.Inlet.CO2(block.Flow2Dir(:,1))); %only makes sense if FCtype=1 and CO2 is a cathode state
    end
    Tags.(block.name).Tpen = T_Elec;
    Tags.(block.name).TcathOut = Flow2Out.T;
    Tags.(block.name).TanodeOut = Flow1Out.T;
    Tags.(block.name).Current = sum(nCurrent);
    Tags.(block.name).StackPower = Tags.(block.name).Current*Tags.(block.name).Voltage*block.Cells/1000; %power in kW
    Tags.(block.name).SinglePassUtilization = Tags.(block.name).Current*block.Cells/(2000*F)/(4*Inlet.Flow1.CH4+Inlet.Flow1.CO + Inlet.Flow1.H2);
    if strcmp(block.CoolingStream,'cathode')
        Tags.(block.name).StackdeltaT = Flow2Out.T-Inlet.Flow2.T;
    elseif strcmp(block.CoolingStream,'anode')
        Tags.(block.name).StackdeltaT = Flow1Out.T-mean(Flow1.Inlet.T(block.Flow1Dir(:,1)));
    end
    Tags.(block.name).PENtemperature = T_Elec;
    Tags.(block.name).PENavgT = sum(T_Elec)/block.nodes;
    Tags.(block.name).MaxPEN = max(T_Elec);
    Tags.(block.name).PENdeltaT = Tags.(block.name).MaxPEN-min(T_Elec);
    Tags.(block.name).dTdX = (T_Elec-T_Elec(block.HTadjacent(:,2)))/(block.L_Cell/block.columns);
    Tags.(block.name).dTdY = (T_Elec-T_Elec(block.HTadjacent(:,4)))/(block.W_Cell/block.rows);
    Tags.(block.name).MaxdTdX = max(abs([Tags.(block.name).dTdX;Tags.(block.name).dTdY;]));
    Tags.(block.name).Efficiency = Tags.(block.name).Power/(802952.15*Inlet.Flow1.CH4 + 305200*Inlet.Flow1.CO + 240424*Inlet.Flow1.H2);
    %% Outlet Ports
    Out.Flow1Out = Flow1Out;
    Out.Flow2Out  = Flow2Out;
    Out.Flow1Pin = P_flow1;
    Out.Flow2Pin = P_flow2;
    Out.MeasureVoltage = Tags.(block.name).Voltage;
    Out.MeasurePower = Tags.(block.name).Power;
    Out.MeasureTflow1 = mean(Y(3*nodes+block.Flow1Dir(:,end)));
    Out.MeasureTflow2 = mean(Y(nodes+block.Flow2Dir(:,end)));
elseif strcmp(string1,'dY')
    Current.CO = Tags.(block.name).I_CO;
    Current.H2 = Tags.(block.name).I_H2;
    [R,FLow1Out] = KineticReformation(block.method,Flow1,P_flow1,block.KineticCoeff1,Current,block);%% Kinetic reaction rates  (WGS is always near equilibrium)
    switch block.Reformer
        case 'internal'
            RefCurrent.H2 = zeros(nodes,1);
            RefCurrent.CO = zeros(nodes,1);
            [Rref,RefOut] = KineticReformation(block.method,Flow3,P_flow1,block.KineticCoeff3,RefCurrent,block);%% Kinetic reaction rates  (WGS is always near equilibrium)
            R.CH4ref = Rref.CH4;
            R.WGSref = Rref.WGS;
    end
    [h,~] = enthalpy(T_Elec,{'H2','H2O','O2','CO','CO2'});
    Power = Tags.(block.name).Voltage*nCurrent/1000; %cell power in kW
    Qgen = Current.H2/(2000*F).*(h.H2+.5*h.O2-h.H2O) + Current.CO/(2000*F).*(h.CO+.5*h.O2-h.CO2) - Power;%kW of heat generated by electrochemistry (per node & per cell)
    Tags.(block.name).Q_gen = sum(Qgen*block.Cells); %kW of heat generated by electrochemistry
    switch block.FCtype%ion transport across membrane (total enthalpy)
        case {'SOFC';'SOEC'}
            Qion = nCurrent/(4000*F).*h.O2; %O2 ion crossing over (kW)
        case {'MCFC';'MCEC'}
            Qion = nCurrent/(4000*F).*h.O2 + nCurrent/(2000*F).*h.CO2;% O2 & CO2 ion crossing over
    end
    %% Q %% Heat transfer & Generation
    switch block.Reformer
        case 'internal'
            QT = block.HTconv*Y(1:6*nodes) + block.HTcond*Y(1:6*nodes) + block.HTrad*(Y(1:6*nodes).^4);
        case 'direct'
            QT = block.HTconv*Y(1:5*nodes) + block.HTcond*Y(1:5*nodes) + block.HTrad*(Y(1:5*nodes).^4);
    end
    
    %energy flows & sepcific heats
    HoutCath = enthalpy(Flow2.Outlet);
    HinCath = enthalpy(Flow2.Inlet);
    HoutAnode = enthalpy(Flow1.Outlet);
    HinAnode = enthalpy(Flow1.Inlet);
    
    %% %% solve for dY in order of states
    dY = 0*Y;
    %%Temperatures
    dY(1:nodes)= QT(1:nodes)./block.tC(1:nodes);  %Ox Sep Plate
    for i=1:1:length(block.Flow2Dir(1,:)) %having the downstream nodes change temperature with the upstream nodes prevents propogation issues when taking larger time steps
        k = block.Flow2Dir(:,i);
        dY(nodes+k)= (QT(nodes+k) + HinCath(k) - HoutCath(k) - Qion(k))./block.tC(nodes+k); %Cathode
        if i>1
            dY(nodes+k) = dY(nodes+k)+dY(nodes+kprev);
        end
        kprev = k;
    end
    dY(1+2*nodes:3*nodes)= (QT(1+2*nodes:3*nodes) + Qgen)./block.tC(2*nodes+1:3*nodes); %Electrolyte Plate
    for i=1:1:length(block.Flow1Dir(1,:))
        k = block.Flow1Dir(:,i);
        dY(3*nodes+k)= (QT(3*nodes+k) + HinAnode(k) - HoutAnode(k) + Qion(k) - Power(k) - Qgen(k))./block.tC(3*nodes+k);  %Anode
        if i>1
            dY(3*nodes+k) = dY(3*nodes+k)+dY(3*nodes+kprev);
        end
        kprev = k;
    end
    dY(1+4*nodes:5*nodes)= QT(1+4*nodes:5*nodes)./block.tC(4*nodes+1:5*nodes);  %Fuel Sep Plate
    
    n =nT;
    %%Cathode Species
    for i = 1:1:length(block.Spec2)
        if strcmp(block.Spec2{i},'CO2') && any(strcmp(block.FCtype,{'MCFC';'MCEC'}))
            dY(n+1:n+nodes)= (Flow2.Inlet.CO2 - Y(n+1:n+nodes) - nCurrent/(2*F*1000))./block.tC(n+1:n+nodes);  %CO2 species concentration with CO2 crossover
        elseif strcmp(block.Spec2{i},'O2') && any(strcmp(block.FCtype,{'SOFC';'SOEC';'MCFC';'MCEC'}))
            dY(n+1:n+nodes)= (Flow2.Inlet.O2 - Y(n+1:n+nodes) - nCurrent/(4*F*1000))./block.tC(n+1:n+nodes);%O2 species concentration with O2 crossover
        else
            dY(n+1:n+nodes)= (Flow2.Inlet.(block.Spec2{i}) - Flow2.Outlet.(block.Spec2{i}))./block.tC(n+1:n+nodes);%all other species concentration
        end
        n = n+nodes;
    end
    %% Anode Species
    for i = 1:1:length(block.Spec1)
        dY(n+1:n+nodes)= (FLow1Out.(block.Spec1{i}) - Flow1.Outlet.(block.Spec1{i}))./block.tC(n+1:n+nodes); %all species concentration
        n = n+nodes; 
    end
    %% Reformer
    switch block.Reformer
        case 'internal'
            Hout3 = enthalpy(Flow3.Outlet);
            Hin3 = enthalpy(Flow3.Inlet);
            for i=1:1:length(block.Flow3Dir(1,:))
                k = block.Flow3Dir(:,i);
                dY(5*nodes+k)= (block.RefSpacing*QT(5*nodes+k) + Hin3(k) - Hout3(k))./block.tC(5*nodes+k);  %Fuel Reformer Channels
                if i>1
                    dY(5*nodes+k) = dY(5*nodes+k)+dY(5*nodes+kprev);
                end
                kprev = k;
            end
            for i = 1:1:length(block.Spec1)
                dY(n+1:n+nodes)= (RefOut.(block.Spec1{i}) - Flow3.Outlet.(block.Spec1{i}))./block.tC(n+1:n+nodes);   %all species concentrations
                n = n+nodes;
            end
    end

    %%Current
    dY(n+1:n+nodes) = (Inlet.NetCurrent - sum(nCurrent) + (Tags.(block.name).nVoltage-Tags.(block.name).Voltage)./Tags.(block.name).ASR.*(block.A_Cell*100^2))./block.tC(end-nodes-1:end-2); n = n+nodes; %error in A/cm^2 * area  
    
    %%Pressure
    Nflow1 = block.Pfactor1*max(0.01,(P_flow1-Inlet.Flow1Pout));%total anode flow out
    dY(n+1) = (NetFlow(Inlet.Flow1)-Nflow1)*Ru*Inlet.Flow1.T/block.tC(n+2);
    if block.ClosedCathode %closed end cathode
        dY(n+2) = 0;%(NetFlow(Inlet.Flow2)-sum(nCurrent)/(4*F*1000))*Ru*Inlet.Flow2.T/block.tC(n+1); %no flow out, so anything not used in electrochemistry adds to pressure
    else
        Nflow2 = block.Pfactor2*max(0.1,(P_flow2-Inlet.Flow2Pout));%total cathode flow out
        dY(n+2) = (NetFlow(Inlet.Flow2)-Nflow2)*Ru*Inlet.Flow2.T/block.tC(n+1);%working with total flow rates so must multiply by nodes & cells
    end
    Out = dY;
end