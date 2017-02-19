function Out = FuelCell(t,Y, Inlet,block,string1)
%This function models a fuel cell  it takes input Y states and returns dY: 
%FC model with many states: Temperatures (oxidizer plate, cathode, electrolyte, anode, fuel plate [reformer]), Cathode species ( [ CO2, H2O], N2, O2) anode species (CH4, CO, CO2, H2, H2O, N2, O2), [Reformer species] [Rate of internal reforming reactions] Current, cathode pressure, anode pressure
% Five (5) inlets: {'PowerDemand','CathodeIn','AnodeIn','CathodePressureOut','AnodePressureOut'}
% Seven (7) outlets: {'CathodeOut','AnodeOut','CathodePressureIn','AnodePressureIn','MeasureCurrent','MeasureTpen','MeasureTcathOut'}
global F Ru Tags
nodes = block.nodes;

%% scale states to physical values
Y = Y.*block.Scale;

%% add species that may not be in inlet
inFields = fieldnames(Inlet.CathodeIn);
for i = 1:1:length(block.CathSpec)
    if ~ismember(block.CathSpec{i},inFields)
        Inlet.CathodeIn.(block.CathSpec{i}) = 0;
    end
end
inFields = fieldnames(Inlet.AnodeIn);
for i = 1:1:length(block.AnSpec)
    if ~ismember(block.AnSpec{i},inFields)
        Inlet.AnodeIn.(block.AnSpec{i}) = 0;
    end
end
%% seperate out temperatures
T_Cath_plate = Y(1:nodes);
Cathode.Outlet.T = Y(nodes+1:2*nodes);
T_Elec = Y(2*nodes+1:3*nodes);
Anode.Outlet.T = Y(3*nodes+1:4*nodes); 
T_Anode_plate = Y(4*nodes+1:5*nodes); 
switch block.Reformer
    case 'internal'
        nT = 6*nodes; % # of temperature states
        Reformer.Outlet.T = Y(5*nodes+1:6*nodes); %only gets used if internal reformer exists, otherwise these values are actually QT1
    case 'adiabatic'
        nT = 5*nodes+1; % # of temperature states
        Reformer.Outlet.T = Y(5*nodes+1);
    case 'direct'
        nT = 5*nodes; % # of temperature states
end
n = nT;

%Current
nCurrent = Y(end-nodes-1:end-2);
%Pressure
Pcath = Y(end-1); %pressure
Panode = Y(end); %pressure

%% Cathode
for i = 1:1:length(block.CathSpec)
    Cathode.Outlet.(block.CathSpec{i}) = Y(n+1:n+nodes); n = n+nodes;
end
for j = 1:1:length(block.AirFlowDir(1,:));
    if j==1%first column recieves fresh inlet
        k = block.AirFlowDir(:,1);
        Cathode.Inlet.T(k,1) = Inlet.CathodeIn.T; 
        for i = 1:1:length(block.CathSpec)
            Cathode.Inlet.(block.CathSpec{i})(k,1) = Inlet.CathodeIn.(block.CathSpec{i})/block.Cells/length(k); 
        end
    else%subsequent columns recieve outlet of previous column
        k2 = block.AirFlowDir(:,j);
        Cathode.Inlet.T(k2,1) = Cathode.Outlet.T(k);
        for i = 1:1:length(block.CathSpec)
            Cathode.Inlet.(block.CathSpec{i})(k2,1) = Cathode.Outlet.(block.CathSpec{i})(k); 
        end
        k = k2;
    end
end
if block.ClosedCathode %closed end cathode
    Cathode.Outlet.O2 = zeros(nodes,1);
    Cathode.Inlet.O2 = zeros(nodes,1);
    k = block.AirFlowDir(:,1);
    Cathode.Inlet.O2(k,1) = Inlet.CathodeIn.O2/block.Cells/length(k); 
    if strcmp(block.FCtype,'MCFC')
        Cathode.Outlet.CO2 = zeros(nodes,1);
        Cathode.Inlet.CO2 = zeros(nodes,1);
        Cathode.Inlet.CO2(k,1) = Inlet.CathodeIn.CO2/block.Cells/length(k); 
    end
    c = length(block.AirFlowDir(1,:));% # of columns
    for j = 1:1:c
        Cathode.Outlet.O2(k,1) = Cathode.Inlet.O2(k,1) - nCurrent(k)/(4*F*1000); 
        if strcmp(block.FCtype,'MCFC')
            Cathode.Outlet.CO2(k,1) = Cathode.Inlet.CO2(k,1) - nCurrent(k)/(2*F*1000); 
        end
        if j<c %subsequent columns recieve outlet of previous column
            k2 = block.AirFlowDir(:,j+1);
            Cathode.Inlet.O2(k2,1) = Cathode.Outlet.O2(k,1); 
            if strcmp(block.FCtype,'MCFC')
                Cathode.Inlet.CO2(k2,1) = Cathode.Outlet.CO2(k,1); 
            end
            k = k2;
        end
    end
end


%% Anode 
for i = 1:1:length(block.AnSpec)
    Anode.Outlet.(block.AnSpec{i}) = Y(n+1:n+nodes); n = n+nodes;
end
%% Internal reformer
switch block.Reformer
    case 'internal'
        for i = 1:1:length(block.AnSpec)
            Reformer.Outlet.(block.AnSpec{i}) = Y(n+1:n+nodes); n = n+nodes;
        end
    case 'adiabatic'
        for i = 1:1:length(block.AnSpec)
            Reformer.Outlet.(block.AnSpec{i}) = Y(n+1); n = n+1;
        end
end

AnodeOut.T  = mean(Anode.Outlet.T(block.FuelFlowDir(:,end))); %temperature 
for i = 1:1:length(block.AnSpec)
    AnodeOut.(block.AnSpec{i}) = max(0,sum(Anode.Outlet.(block.AnSpec{i})(block.FuelFlowDir(:,end)))*block.Cells);%avoid sending negative outlets
end
CathodeOut.T  = mean(Cathode.Outlet.T(block.AirFlowDir(:,end))); %temperature 
for i = 1:1:length(block.CathSpec)
    CathodeOut.(block.CathSpec{i}) = max(0,sum(Cathode.Outlet.(block.CathSpec{i})(block.AirFlowDir(:,end)))*block.Cells);%avoid sending negative outlets
end
%% Reformer
switch block.Reformer
    case 'internal'
        for j = 1:1:length(block.ReformFlowDir(1,:))
            if j==1
                k = block.ReformFlowDir(:,1);
                Reformer.Inlet.T(k,1) = Inlet.AnodeIn.T;
                for i = 1:1:length(block.AnSpec)
                    Reformer.Inlet.(block.AnSpec{i})(k,1) = Inlet.AnodeIn.(block.AnSpec{i})/block.Cells/length(k)*block.RefSpacing;
                end
            else
                k2 = block.ReformFlowDir(:,j);
                Reformer.Inlet.T(k2,1) = Reformer.Outlet.T(k);
                for i = 1:1:length(block.AnSpec)
                    Reformer.Inlet.(block.AnSpec{i})(k2,1) = Reformer.Outlet.(block.AnSpec{i})(k);
                end
                k = k2;
            end
        end
    case 'adiabatic'
        Reformer.Inlet = Inlet.AnodeIn;
end
%% Anode
for j = 1:1:length(block.FuelFlowDir(1,:))
    k2 =block.FuelFlowDir(:,j);
    if j==1 % first column of fuel flow direction
        switch block.Reformer
            case 'internal'
                Anode.Inlet.T(k2,1) =  Reformer.Outlet.T(k);
            case 'adiabatic'
                Anode.Inlet.T(k2,1) = Reformer.Outlet.T;
            case 'direct'
                Anode.Inlet.T(k2,1) = Inlet.AnodeIn.T;
        end
        for i = 1:1:length(block.AnSpec)
            switch block.Reformer
                case 'internal'
                    Anode.Inlet.(block.AnSpec{i})(k2,1) = Reformer.Outlet.(block.AnSpec{i})(k)/block.RefSpacing;%Species flows coming into anode
                case 'adiabatic'
                    Anode.Inlet.(block.AnSpec{i})(k2,1) = Reformer.Outlet.(block.AnSpec{i})/block.Cells/length(k2);
                case 'direct'
                    Anode.Inlet.(block.AnSpec{i})(k2,1) = Inlet.AnodeIn.(block.AnSpec{i})/block.Cells/length(k2);
            end
        end
    else
        Anode.Inlet.T(k2,1) = Anode.Outlet.T(k);
        for i = 1:1:length(block.AnSpec)
            Anode.Inlet.(block.AnSpec{i})(k2,1) = Anode.Outlet.(block.AnSpec{i})(k);
        end
    end
    k = k2;
end

if strcmp(string1,'Outlet')
    %% Nernst & Losses
    n_an_in = NetFlow(Anode.Inlet);
    n_an_out = NetFlow(Anode.Outlet);
    n_cath_in = NetFlow(Cathode.Inlet);
    n_cath_out = NetFlow(Cathode.Outlet);
    switch block.FCtype
        case 'SOFC'
            if block.ClosedCathode
                AvgX.O2 = ones(block.nodes,1);
            else
                AvgX.O2 = (Cathode.Outlet.O2+Cathode.Inlet.O2)./(n_cath_in+n_cath_out);
            end
        case 'MCFC'
            if block.ClosedCathode
                AvgX.O2 = ones(block.nodes,1)/3;
                AvgX.CO2c = 2*ones(block.nodes,1)/3;
            else
                AvgX.O2 = (Cathode.Outlet.O2+Cathode.Inlet.O2)./(n_cath_in+n_cath_out);
                AvgX.CO2c = (Cathode.Outlet.CO2+Cathode.Inlet.CO2)./(n_cath_in+n_cath_out);
            end
            AvgX.CO2a = (Anode.Outlet.CO2+Anode.Inlet.CO2)./(n_an_in+n_an_out);
    end

    AvgX.H2 = (Anode.Outlet.H2+Anode.Inlet.H2)./(n_an_in+n_an_out);
    AvgX.H2O = (Anode.Outlet.H2O+Anode.Inlet.H2O)./(n_an_in+n_an_out);
    
    k = block.FuelFlowDir(:,1);
    if min(Anode.Inlet.H2(k))==0 %gives some reformed methane as anode inlet
        R.CH4 = Anode.Inlet.CH4 - Anode.Outlet.CH4; %Rate of methane reforming R.CH4
        K_WGS = exp(4189.8./Anode.Outlet.T -3.8242);% Water gas shift equilibrium constant
        CO_eq = Anode.Outlet.CO2.*Anode.Outlet.H2./(K_WGS.*Anode.Outlet.H2O);
        R.WGS = (Anode.Inlet.CO+R.CH4)-CO_eq; %inlet CO + CO from reforming - outlet CO
        AvgX.H2 = (Anode.Outlet.H2 + 0.5*(3*R.CH4(k)+R.WGS(k)))./(n_an_in(k)+n_an_out(k));
        AvgX.H2O = (Anode.Outlet.H2O - 0.5*(R.CH4(k) + R.WGS(k))+ Anode.Inlet.H2O(k))./(n_an_in(k)+n_an_out(k));
    end
    FuelCellNernst(nCurrent,T_Elec,Pcath,AvgX,block);
    Voltage =  sum(Tags.(block.name).nVoltage.*(nCurrent/sum(nCurrent)));
    Tags.(block.name).Voltage = Voltage;
    %% Tags
    H2_in = sum(4*Anode.Inlet.CH4(block.FuelFlowDir(:,1))+ Anode.Inlet.CO(block.FuelFlowDir(:,1)) + Anode.Inlet.H2(block.FuelFlowDir(:,1)));
    H2_out = sum(4*Anode.Outlet.CH4(block.FuelFlowDir(:,end))+ Anode.Outlet.CO(block.FuelFlowDir(:,end)) + Anode.Outlet.H2(block.FuelFlowDir(:,end)));
    Tags.(block.name).H2utilization = (H2_in - H2_out)./H2_in;
    Tags.(block.name).O2utilization = sum(Cathode.Inlet.O2(block.AirFlowDir(:,1)) - Cathode.Outlet.O2(block.AirFlowDir(:,end)))/sum(Cathode.Inlet.O2(block.AirFlowDir(:,1)));
    if strcmp(block.FCtype,'MCFC')
        Tags.(block.name).CO2utilization = sum(Cathode.Inlet.CO2(block.AirFlowDir(:,1)) - Cathode.Outlet.CO2(block.AirFlowDir(:,end)))/sum(Cathode.Inlet.CO2(block.AirFlowDir(:,1))); %only makes sense if FCtype=1 and CO2 is a cathode state
    end
    Tags.(block.name).Tpen = T_Elec;
    Tags.(block.name).TcathOut = CathodeOut.T;
    Tags.(block.name).TanodeOut = AnodeOut.T;
    Tags.(block.name).Current = sum(nCurrent);
    if strcmp(block.CoolingStream,'cathode')
        Tags.(block.name).StackdeltaT = CathodeOut.T-Inlet.CathodeIn.T;
    elseif strcmp(block.CoolingStream,'anode')
        Tags.(block.name).StackdeltaT = AnodeOut.T-mean(Anode.Inlet.T(block.FuelFlowDir(:,1)));
    end
    Tags.(block.name).PENavgT = sum(T_Elec)/block.nodes;
    Tags.(block.name).MaxPEN = max(T_Elec);
    Tags.(block.name).PENdeltaT = Tags.(block.name).MaxPEN-min(T_Elec);
    Tags.(block.name).dTdX = (T_Elec-T_Elec(block.HTadjacent(:,2)))/(block.L_Cell/block.columns);
    Tags.(block.name).dTdY = (T_Elec-T_Elec(block.HTadjacent(:,4)))/(block.W_Cell/block.rows);
    Tags.(block.name).MaxdTdX = max(abs([Tags.(block.name).dTdX;Tags.(block.name).dTdY;]));
    Tags.(block.name).AnCH4last = Anode.Outlet.CH4(block.FuelFlowDir(:,end));
    %% Outlet Ports
    Out.CathodeOut  = CathodeOut;
    Out.AnodeOut = AnodeOut;
    Out.CathodePressureIn = Pcath;
    Out.AnodePressureIn = Panode;
    Out.MeasureVoltage = Voltage;
    Out.MeasurePower = Tags.(block.name).Power;
    Out.MeasureTpen = Y(2*nodes+1:3*nodes);
    Out.MeasureTcathOut = Y(nodes+block.AirFlowDir(:,end));
    Out.MeasureTanodeOut = Y(3*nodes+block.FuelFlowDir(:,end));
elseif strcmp(string1,'dY')  
    Voltage = Tags.(block.name).Voltage;
    [R,AnOut] = KineticReformation(block.method,Anode,Panode,block.KineticCoeff1,nCurrent,block);%% Kinetic reaction rates  (WGS is always near equilibrium)
    switch block.Reformer
        case 'internal'
            [Rref,RefOut] = KineticReformation(block.method,Reformer,Panode,block.KineticCoeff3,zeros(nodes,1),block);%% Kinetic reaction rates  (WGS is always near equilibrium)
            R.CH4ref = Rref.CH4;
            R.WGSref = Rref.WGS;
%         case 'adiabatic'
%             [Rref,RefOut] = KineticReformation('equilibrium',Reformer,Panode,0,0,block);%% Kinetic reaction rates  (WGS is always near equilibrium)
%             R.CH4ref = Rref.CH4;
%             R.WGSref = Rref.WGS;
    end
    [h,~] = enthalpy(T_Elec,{'H2','H2O','O2','CO2'});
    Power = Voltage*nCurrent/1000; %cell power in kW
    Qgen = nCurrent/(2000*F).*(h.H2+.5*h.O2-h.H2O)-Power;%kW of heat generated by electrochemistry (per node & per cell)
    switch block.FCtype%ion transport across membrane (total enthalpy)
        case 'SOFC'
            Qion = nCurrent/(4000*F).*h.O2; %O2 ion crossing over (kW)
        case 'MCFC'
            Qion = nCurrent/(4000*F).*h.O2 + nCurrent/(2000*F).*h.CO2;% O2 & CO2 ion crossing over
    end
    %% Q %% Heat transfer & Generation
    switch block.Reformer
        case 'internal'
            QT = block.HTmatrix*Y(1:6*nodes);
        case {'adiabatic';'direct';'none';'external';}
            QT = block.HTmatrix*Y(1:5*nodes);
    end
%     QT = QT + block.RTmatrix*Y(1:length(QT)).^4;
    Tags.(block.name).Q_gen = sum(Qgen*block.Cells); %kW of heat generated by electrochemistry
    %energy flows & sepcific heats
    HoutCath = enthalpy(Cathode.Outlet);
    HinCath = enthalpy(Cathode.Inlet);
    HoutAnode = enthalpy(Anode.Outlet);
    HinAnode = enthalpy(Anode.Inlet);
    
    %% %% solve for dY in order of states
    dY = 0*Y;
    %%Temperatures
    dY(1:nodes)= QT(1:nodes)./block.tC(1:nodes);  %Ox Sep Plate
    for i=1:1:length(block.AirFlowDir(1,:)) %having the downstream nodes change temperature with the upstream nodes prevents propogation issues when taking larger time steps
        k = block.AirFlowDir(:,i);
        dY(nodes+k)= (QT(nodes+k) + HinCath(k) - HoutCath(k) - Qion(k))./block.tC(nodes+k); %Cathode
        if i>1
            dY(nodes+k) = dY(nodes+k)+dY(nodes+kprev);
        end
        kprev = k;
    end
    dY(1+2*nodes:3*nodes)= (QT(1+2*nodes:3*nodes) + Qgen)./block.tC(2*nodes+1:3*nodes); %Electrolyte Plate
    for i=1:1:length(block.FuelFlowDir(1,:))
        k = block.FuelFlowDir(:,i);
        dY(3*nodes+k)= (QT(3*nodes+k) + HinAnode(k) - HoutAnode(k) + Qion(k) - Power(k) - Qgen(k))./block.tC(3*nodes+k);  %Anode
        if i>1
            dY(3*nodes+k) = dY(3*nodes+k)+dY(3*nodes+kprev);
        end
        kprev = k;
    end
    dY(1+4*nodes:5*nodes)= QT(1+4*nodes:5*nodes)./block.tC(4*nodes+1:5*nodes);  %Fuel Sep Plate
    
    
    n =nT;
    %%Cathode Species
    for i = 1:1:length(block.CathSpec)
        if strcmp(block.CathSpec{i},'CO2') && strcmp(block.FCtype,'MCFC')
            dY(n+1:n+nodes)= (Cathode.Inlet.CO2 - Cathode.Outlet.CO2 - nCurrent/(2*F*1000))./block.tC(n+1:n+nodes);  %CO2 species concentration with CO2 crossover
        elseif strcmp(block.CathSpec{i},'O2') && (strcmp(block.FCtype,'SOFC') || strcmp(block.FCtype,'MCFC'))
            dY(n+1:n+nodes)= (Cathode.Inlet.O2 - Cathode.Outlet.O2 - nCurrent/(4*F*1000))./block.tC(n+1:n+nodes);%O2 species concentration with O2 crossover
        else
            dY(n+1:n+nodes)= (Cathode.Inlet.(block.CathSpec{i}) - Cathode.Outlet.(block.CathSpec{i}))./block.tC(n+1:n+nodes);%all other species concentration
        end
        n = n+nodes;
    end
    %% Anode Species
    for i = 1:1:length(block.AnSpec)
        dY(n+1:n+nodes)= (AnOut.(block.AnSpec{i}) - Anode.Outlet.(block.AnSpec{i}))./block.tC(n+1:n+nodes); %all species concentration
        n = n+nodes; 
    end
   
    %% Reformer
    switch block.Reformer
        case 'internal'
            HoutReform = enthalpy(Reformer.Outlet);
            HinReform = enthalpy(Reformer.Inlet);
            for i=1:1:length(block.ReformFlowDir(1,:))
                k = block.ReformFlowDir(:,i);
                dY(5*nodes+k)= (block.RefSpacing*QT(5*nodes+k) + HinReform(k) - HoutReform(k))./block.tC(5*nodes+k);  %Fuel Reformer Channels
                if i>1
                    dY(5*nodes+k) = dY(5*nodes+k)+dY(5*nodes+kprev);
                end
                kprev = k;
            end
             for i = 1:1:length(block.AnSpec)
                dY(n+1:n+nodes)= (RefOut.(block.AnSpec{i}) - Reformer.Outlet.(block.AnSpec{i}))./block.tC(n+1:n+nodes);   %all species concentrations
                n = n+nodes;
             end
%         case 'adiabatic'
%             dY(5*nodes+1)= (enthalpy(Reformer.Inlet)-enthalpy(Reformer.Outlet))/block.tC(5*nodes+1);  %Fuel Reformer
%              for i = 1:1:length(block.AnSpec)
%                  dY(n+1)= (RefOut.(block.AnSpec{i}) - Reformer.Outlet.(block.AnSpec{i}))./block.tC(n+1);  %all species concentrations
%                  n = n+1; 
%              end
    end

    %%Current
    dY(n+1:n+nodes) = (Inlet.NetCurrent - sum(nCurrent) + (Tags.(block.name).nVoltage-Voltage)./Tags.(block.name).ASR.*(block.A_Cell*100^2))./block.tC(end-nodes-1:end-2); n = n+nodes; %error in A/cm^2 * area  
    Power = Voltage*Inlet.NetCurrent*block.Cells/1000;
    Efficiency = Power/(802952.15*Inlet.AnodeIn.CH4 + 305200*Inlet.AnodeIn.CO + 240424*Inlet.AnodeIn.H2);
    Tags.(block.name).Efficiency = Efficiency;

    %% Pressure
    Nanode = block.PfactorAnode*max(0.01,(Panode-Inlet.AnodePressureOut));%total anode flow out
    if block.ClosedCathode %closed end cathode
        dY(n+1) = 0;%(NetFlow(Inlet.CathodeIn)-sum(nCurrent)/(4*F*1000))*Ru*Inlet.CathodeIn.T/block.tC(n+1); %no flow out, so anything not used in electrochemistry adds to pressure
    else
        Ncath = block.PfactorCath*max(0.1,(Pcath-Inlet.CathodePressureOut));%total cathode flow out
        dY(n+1) = (NetFlow(Inlet.CathodeIn)-Ncath)*Ru*Inlet.CathodeIn.T/block.tC(n+1);%working with total flow rates so must multiply by nodes & cells
    end
    dY(n+2) = (NetFlow(Inlet.AnodeIn)-Nanode)*Ru*Inlet.AnodeIn.T/block.tC(n+2);

    dY = dY./block.Scale; 
    Out = dY;
end

function [R,Out] = KineticReformation(method,Flow,Pressure,KineticCoeff,Current,block)%find new reforming reaction rates
global Ru F
spec = block.AnSpec;
Type = block.FCtype;
K_WGS = exp(4189.8./Flow.Outlet.T -3.8242);% Water gas shift equilibrium constant
nout = NetFlow(Flow.Outlet);
% nin = NetFlow(Flow.Inlet);
% X_CH4 = (Flow.Inlet.CH4+Flow.Outlet.CH4)./(nin+nout)*Pressure*1000; %partial pressures in Pa
% X_H2O = (Flow.Inlet.H2O+Flow.Outlet.H2O)./(nin+nout)*Pressure*1000; %partial pressures in Pa
    X_CH4 = Flow.Outlet.CH4./nout*Pressure*1000; %partial pressures in Pa
    X_H2O = Flow.Outlet.H2O./nout*Pressure*1000; %partial pressures in Pa
if strcmp(method,'Achenbach')
    R.CH4 = KineticCoeff*(X_CH4.*exp(-8.2e4./(Ru*Flow.Outlet.T)));
elseif strcmp(method,'Leinfelder')
    R.CH4 = KineticCoeff*(30.8e10*X_CH4.*X_H2O.*exp(-2.05e5./(Ru*Flow.Outlet.T)));
elseif strcmp(method,'Drescher')
    R.CH4 = KineticCoeff*(288.52*X_CH4.*X_H2O.*exp(-1.1e4./(Ru*Flow.Outlet.T))./(1+16*X_CH4+0.143*X_H2O.*exp(3.9e4./(Ru*Flow.Outlet.T))));
elseif strcmp(method,'equilibrium')
    a = 4352.2./Flow.Outlet.T - 3.99;
    K_WGS = block.scaleK_WGS.*exp(a);% Water gas shift equilibrium constant
    K_CH4 = block.scaleK_CH4.*2459000.*exp(-6.187*a);
    CH4_eq = (Flow.Outlet.H2.^3.*Flow.Outlet.CO)./(K_CH4.*Flow.Outlet.H2O).*(Pressure./NetFlow(Flow.Outlet)).^2;
    R.CH4 = Flow.Inlet.CH4-CH4_eq;
end
%confirm we don't violate anything casuing negative species
R.CH4 = min([R.CH4,Flow.Inlet.CH4,Flow.Inlet.H2O],[],2);
R.CH4 = max([R.CH4,-Flow.Inlet.CO,-Flow.Inlet.H2/3],[],2);

CO_eq = Flow.Outlet.CO2.*Flow.Outlet.H2./(K_WGS.*Flow.Outlet.H2O);
R.WGS = (Flow.Inlet.CO+R.CH4)-CO_eq; %inlet CO + CO from reforming - outlet CO
R.WGS = min([R.WGS,Flow.Inlet.CO+R.CH4,Flow.Inlet.H2O-R.CH4],[],2);
R.WGS = max([R.WGS,-Flow.Inlet.CO2,-(Flow.Inlet.H2+3*R.CH4)],[],2);
%identify the target outflow that dY will converge to
for i = 1:1:length(spec)
    if strcmp(spec{i},'CO2') 
        Out.CO2 = (Flow.Inlet.CO2 + R.WGS); 
        if strcmp(Type,'MCFC')
            Out.CO2 = Out.CO2 + Current/(2*F*1000);
        end
    elseif strcmp(spec{i},'H2') 
        Out.H2 = (Flow.Inlet.H2 + 3*R.CH4 + R.WGS - Current/(2*F*1000));
    elseif strcmp(spec{i},'H2O') 
        Out.H2O = (Flow.Inlet.H2O - R.CH4 - R.WGS + Current/(2*F*1000)); 
    elseif strcmp(spec{i},'CH4') 
        Out.CH4 = (Flow.Inlet.CH4 - R.CH4); 
    elseif strcmp(spec{i},'CO') 
        Out.CO = (Flow.Inlet.CO + R.CH4 - R.WGS);  
    else
        Out.(spec{i}) = Flow.Inlet.(spec{i});  
    end
end 