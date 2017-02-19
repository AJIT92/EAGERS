function Out = Electrolyzer(t,Y, Inlet,block,string1)
%This function models an electorlyzer it takes input Y states and returns dY: 
%Electrolyzer model with many states: Temperatures (oxidizer plate, cathode, electrolyte, anode, fuel plate [reformer]), Cathode species ( [ CO2, H2O], N2, O2) anode species (CH4, CO, CO2, H2, H2O, N2, O2), [Reformer species] [Rate of internal reforming reactions] Current, cathode pressure, anode pressure
% Five (5) inlets: {'PowerDemand','CathodeIn','AnodeIn','CathodePressureOut','AnodePressureOut'}
% Seven (7) outlets: {'CathodeOut','AnodeOut','CathodePressureIn','AnodePressureIn','MeasureVoltage','MeasureTpen','MeasureTcathOut'}
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
    case 'methanator'
        nT = 6*nodes; % # of temperature states
        Methanator.Outlet.T = Y(5*nodes+1:6*nodes); %only gets used if internal Methanator exists, otherwise these values are actually QT1
    case 'none'
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
for j = 1:1:length(block.CathodeFlowDir(1,:));
    if j==1%first column recieves fresh inlet
        k = block.CathodeFlowDir(:,1);
        Cathode.Inlet.T(k,1) = Inlet.CathodeIn.T; 
        for i = 1:1:length(block.CathSpec)
            Cathode.Inlet.(block.CathSpec{i})(k,1) = Inlet.CathodeIn.(block.CathSpec{i})/block.Cells/length(k); 
        end
    else%subsequent columns recieve outlet of previous column
        k2 = block.CathodeFlowDir(:,j);
        Cathode.Inlet.T(k2,1) = Cathode.Outlet.T(k);
        for i = 1:1:length(block.CathSpec)
            Cathode.Inlet.(block.CathSpec{i})(k2,1) = Cathode.Outlet.(block.CathSpec{i})(k); 
        end
        k = k2;
    end
end
CathodeOut.T  = mean(Cathode.Outlet.T(block.CathodeFlowDir(:,end))); %temperature 
for i = 1:1:length(block.CathSpec)
    CathodeOut.(block.CathSpec{i}) = max(0,sum(Cathode.Outlet.(block.CathSpec{i})(block.CathodeFlowDir(:,end)))*block.Cells);%avoid sending negative outlets
end

%% Anode 
for i = 1:1:length(block.AnSpec)
    Anode.Outlet.(block.AnSpec{i}) = Y(n+1:n+nodes); n = n+nodes;
end
for j = 1:1:length(block.AnodeFlowDir(1,:))
    k2 =block.AnodeFlowDir(:,j);
    if j==1 % first column of fuel flow direction
        Anode.Inlet.T(k2,1) = Inlet.AnodeIn.T;
        for i = 1:1:length(block.AnSpec)
            Anode.Inlet.(block.AnSpec{i})(k2,1) = Inlet.AnodeIn.(block.AnSpec{i})/block.Cells/length(k2);
        end
    else
        Anode.Inlet.T(k2,1) = Anode.Outlet.T(k);
        for i = 1:1:length(block.AnSpec)
            Anode.Inlet.(block.AnSpec{i})(k2,1) = Anode.Outlet.(block.AnSpec{i})(k);
        end
    end
    k = k2;
end
AnodeOut.T  = mean(Anode.Outlet.T(block.AnodeFlowDir(:,end))); %temperature 
for i = 1:1:length(block.AnSpec)
    AnodeOut.(block.AnSpec{i}) = max(0,sum(Anode.Outlet.(block.AnSpec{i})(block.AnodeFlowDir(:,end)))*block.Cells);%avoid sending negative outlets
end

%% Methanator
switch block.Reformer
    case 'methanator' %secondary inlet introduces CO2 stream
        for i = 1:1:length(block.AnSpec)
            Methanator.Outlet.(block.AnSpec{i}) = Y(n+1:n+nodes); n = n+nodes;
        end
        for j = 1:1:length(block.MethanatorFlowDir(1,:))
            if j==1
                k = block.MethanatorFlowDir(:,1);
                Methanator.Inlet.T(k,1) = Inlet.CarbonDioxide.T;
                for i = 1:1:length(block.AnSpec)
                    Methanator.Inlet.(block.AnSpec{i})(k,1) = Anode.Outlet.(block.AnSpec{i})(block.AnodeFlowDir(:,1))*block.MethSpacing;
                    if isfield(Inlet.CarbonDioxide,block.AnSpec{i})
                        Methanator.Inlet.(block.AnSpec{i})(k,1) = Methanator.Inlet.(block.AnSpec{i})(k,1) +Inlet.CarbonDioxide.(block.AnSpec{i})/block.Cells/length(k)*block.MethSpacing;
                    end
                end
            else
                k2 = block.MethanatorFlowDir(:,j);
                Methanator.Inlet.T(k2,1) = Methanator.Outlet.T(k);
                for i = 1:1:length(block.AnSpec)
                    Methanator.Inlet.(block.AnSpec{i})(k2,1) = Methanator.Outlet.(block.AnSpec{i})(k);
                end
                k = k2;
            end
        end
end
if strcmp(string1,'Outlet')
    %% Nernst & Losses
    n_an_in = NetFlow(Anode.Inlet);
    n_an_out = NetFlow(Anode.Outlet);
    n_cath_in = NetFlow(Cathode.Inlet);
    n_cath_out = NetFlow(Cathode.Outlet);
    AvgX.O2 = (Anode.Outlet.O2+Anode.Inlet.O2)./(n_an_in+n_an_out);
    AvgX.H2 = (Cathode.Outlet.H2+Cathode.Inlet.H2)./(n_cath_in+n_cath_out);
    AvgX.H2O = (Cathode.Outlet.H2O+Cathode.Inlet.H2O)./(n_cath_in+n_cath_out);

    FuelCellNernst(nCurrent,T_Elec,Panode,AvgX,block);
    Voltage =  sum(Tags.(block.name).nVoltage.*(nCurrent/sum(nCurrent)));
    Tags.(block.name).Voltage = Voltage;
    %% Tags
    H2O_in = Inlet.CathodeIn.H2O;
    H2O_out = CathodeOut.H2O;
    Tags.(block.name).H2Outilization = (H2O_in - H2O_out)./H2O_in;
    Tags.(block.name).Tpen = T_Elec;
    Tags.(block.name).TcathOut = CathodeOut.T;
    Tags.(block.name).TanodeOut = AnodeOut.T;
    Tags.(block.name).Current = sum(abs(nCurrent));
    Tags.(block.name).StackdeltaT = AnodeOut.T-Inlet.AnodeIn.T;
    Tags.(block.name).PENavgT = sum(T_Elec)/block.nodes;
    Tags.(block.name).MaxPEN = max(T_Elec);
    Tags.(block.name).PENdeltaT = Tags.(block.name).MaxPEN-min(T_Elec);
    Tags.(block.name).dTdX = (T_Elec-T_Elec(block.HTadjacent(:,2)))/(block.L_Cell/block.columns);
    Tags.(block.name).dTdY = (T_Elec-T_Elec(block.HTadjacent(:,4)))/(block.W_Cell/block.rows);
    Tags.(block.name).MaxdTdX = max(abs([Tags.(block.name).dTdX;Tags.(block.name).dTdY;]));
    %% Outlet Ports
    Out.CathodeOut  = CathodeOut;
    Out.AnodeOut = AnodeOut;
    Out.CathodePressureIn = Pcath;
    Out.AnodePressureIn = Panode;
    Out.MeasureVoltage = Voltage;
    Out.MeasurePower = Tags.(block.name).Power;
    Out.MeasureTpen = Y(2*nodes+1:3*nodes);
    Out.MeasureTcathOut = Y(nodes+block.CathodeFlowDir(:,end));
    Out.MeasureTanodeOut = Y(3*nodes+block.AnodeFlowDir(:,end));
elseif strcmp(string1,'dY')  
    Voltage =  Tags.(block.name).Voltage;
    switch block.Reformer
        case 'methanator'
            [~,RefOut] = KineticReformation(block.method,Methanator,Panode,block.KineticCoeff3,block.AnSpec);%% Kinetic reaction rates  (WGS is always near equilibrium)
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
            QT = block.HTmatrix*Y(1:6*nodes);
        case {'none'}
            QT = block.HTmatrix*Y(1:5*nodes);
    end
    Tags.(block.name).Q_gen = sum(Qgen*block.Cells); %kW of heat generated by electrochemistry
    %energy flows & sepcific heats
    HoutCath = enthalpy(Cathode.Outlet);
    HinCath = enthalpy(Cathode.Inlet);
    HoutAnode = enthalpy(Anode.Outlet);
    HinAnode = enthalpy(Anode.Inlet);

    %% %% solve for dY in order of states
    dY = 0*Y;
    %%Temperatures
    dY(1:nodes)= QT(1:nodes)./block.tC(1:nodes);  %Cathode Plate
    for i=1:1:length(block.CathodeFlowDir(1,:)) %having the downstream nodes change temperature with the upstream nodes prevents propogation issues when taking larger time steps
        k = block.CathodeFlowDir(:,i);
        dY(nodes+k)= (QT(nodes+k) + HinCath(k) - HoutCath(k)  - Power(k) - Qgen(k) + Qion(k))./block.tC(nodes+k); %Cathode
        if i>1
            dY(nodes+k) = dY(nodes+k)+dY(nodes+kprev);
        end
        kprev = k;
    end
    dY(1+2*nodes:3*nodes)= (QT(1+2*nodes:3*nodes) + Qgen)./block.tC(2*nodes+1:3*nodes); %Electrolyte Plate
    for i=1:1:length(block.AnodeFlowDir(1,:))
        k = block.AnodeFlowDir(:,i);
        dY(3*nodes+k)= (QT(3*nodes+k) + HinAnode(k) - HoutAnode(k) - Qion(k))./block.tC(3*nodes+k);  %Anode
        if i>1
            dY(3*nodes+k) = dY(3*nodes+k)+dY(3*nodes+kprev);
        end
        kprev = k;
    end
    dY(1+4*nodes:5*nodes)= QT(1+4*nodes:5*nodes)./block.tC(4*nodes+1:5*nodes);  %Anode Plate

    n =nT;
    %% Cathode Species
    for i = 1:1:length(block.CathSpec)
        if strcmp(block.CathSpec{i},'H2O')
            dY(n+1:n+nodes)= (Cathode.Inlet.H2O - Cathode.Outlet.H2O + nCurrent/(2*F*1000))./block.tC(n+1:n+nodes);  %H2O species concentration with CO2 crossover
        elseif strcmp(block.CathSpec{i},'H2') 
            dY(n+1:n+nodes)= (Cathode.Inlet.H2 - Cathode.Outlet.H2 - nCurrent/(2*F*1000))./block.tC(n+1:n+nodes);%H2 species concentration with O2 crossover
        else
            dY(n+1:n+nodes)= (Cathode.Inlet.(block.CathSpec{i}) - Cathode.Outlet.(block.CathSpec{i}))./block.tC(n+1:n+nodes);%all other species concentration
        end
        n = n+nodes;
    end

    %% Anode Species
    for i = 1:1:length(block.AnSpec)
        if strcmp(block.AnSpec{i},'O2')
            dY(n+1:n+nodes)= (Anode.Inlet.O2 - Anode.Outlet.O2 - nCurrent/(4*F*1000))./block.tC(n+1:n+nodes);%O2 species concentration with O2 crossover
        else
            dY(n+1:n+nodes)= (Anode.Inlet.(block.AnSpec{i}) - Anode.Outlet.(block.AnSpec{i}))./block.tC(n+1:n+nodes); %all species concentration
        end
        n = n+nodes; 
    end
   
    %% Methanator
    switch block.Reformer
        case 'methanator'
            [~,HoutReform] = enthalpy(Methanator.Outlet);
            [~,HinReform] = enthalpy(Methanator.Inlet);
            for i=1:1:length(block.ReformFlowDir(1,:))
                k = block.ReformFlowDir(:,i);
                dY(5*nodes+k)= (block.RefSpacing*QT(5*nodes+k) + HinReform(k) - HoutReform(k))./block.tC(5*nodes+k);  %Fuel Methanator Channels
                if i>1
                    dY(5*nodes+k) = dY(5*nodes+k)+dY(5*nodes+kprev);
                end
                kprev = k;
            end
            for i = 1:1:length(block.AnSpec)
                dY(n+1:n+nodes)= (RefOut.(block.AnSpec{i}) - Methanator.Outlet.(block.AnSpec{i}))./block.tC(n+1:n+nodes);   %all species concentrations
                n = n+nodes;
            end
    end
    %%Current % note sign convention of current is negative for electrolyzers !!
    dY(n+1:n+nodes) = (Inlet.NetCurrent - sum(nCurrent) + (Tags.(block.name).nVoltage - Voltage)./Tags.(block.name).ASR.*(block.A_Cell*100^2))./block.tC(end-nodes-1:end-2); n = n+nodes; %error in A/cm^2 * area  
    Power = Voltage*abs(Inlet.NetCurrent)*block.Cells/1000;
    Efficiency = 240424*CathodeOut.H2/Power;
    Tags.(block.name).Efficiency = Efficiency;

    %% Pressure
    Nanode = block.PfactorAnode*max(0.01,(Panode-Inlet.AnodePressureOut));%total anode flow out
    Ncath = block.PfactorCath*max(0.1,(Pcath-Inlet.CathodePressureOut));%total cathode flow out
    dY(n+1) = (NetFlow(Inlet.CathodeIn)-Ncath)*Ru*Inlet.CathodeIn.T/block.tC(n+1);%working with total flow rates so must multiply by nodes & cells
    dY(n+2) = (NetFlow(Inlet.AnodeIn)-Nanode)*Ru*Inlet.AnodeIn.T/block.tC(n+2);

    dY = dY./block.Scale; 
    Out = dY;
end

function [R,Out] = KineticReformation(method,Flow,Pressure,KineticCoeff,spec)%find new reforming reaction rates
global Ru
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
end
% a = 4352.2./Flow.Outlet.T - 3.99;
% K_WGS = block.scaleK_WGS.*exp(a);% Water gas shift equilibrium constant
K_WGS = exp(4189.8./Flow.Outlet.T -3.8242);% Water gas shift equilibrium constant
CO_eq = Flow.Outlet.CO2.*Flow.Outlet.H2./(K_WGS.*Flow.Outlet.H2O);
R.WGS = (Flow.Inlet.CO+R.CH4)-CO_eq; %inlet CO + CO from reforming - outlet CO
%confirm we don't violate anything casuing negative species
R.CH4 = min([R.CH4,Flow.Inlet.CH4,Flow.Inlet.H2O],[],2);
R.CH4 = max([R.CH4,-Flow.Inlet.CO,-Flow.Inlet.H2/3],[],2);
R.WGS = min([R.WGS,Flow.Inlet.CO+R.CH4,Flow.Inlet.H2O-R.CH4],[],2);
R.WGS = max([R.WGS,-Flow.Inlet.CO2,-(Flow.Inlet.H2+3*R.CH4)],[],2);
%identify the target outflow that dY will converge to
for i = 1:1:length(spec)
    if strcmp(spec{i},'CO2') 
        Out.CO2 = Flow.Inlet.CO2 + R.WGS; 
    elseif strcmp(spec{i},'H2') 
        Out.H2 = Flow.Inlet.H2 + 3*R.CH4 + R.WGS;
    elseif strcmp(spec{i},'H2O') 
        Out.H2O = Flow.Inlet.H2O - R.CH4 - R.WGS; 
    elseif strcmp(spec{i},'CH4') 
        Out.CH4 = Flow.Inlet.CH4 - R.CH4; 
    elseif strcmp(spec{i},'CO') 
        Out.CO = Flow.Inlet.CO + R.CH4 - R.WGS;  
    else
        Out.(spec{i}) = Flow.Inlet.(spec{i});  
    end
end 