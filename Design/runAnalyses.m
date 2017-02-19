function Out = runAnalyses(Project,varargin)
MSGra = msgbox('Re-calculating Results');
demand=Project.Building; %the building load structure
if isfield(Project,'Renewable')
    renew = Project.Renewable;
else renew = 1;
end
if ~isempty(varargin)
    user = varargin{1};
    scale = varargin{2};
else user = 'Commercial';
    scale = 1;
end

%% Dispatch System
Out=dispatchSystem(Project.System,demand,Project.Utilities.Grid,renew,Project.Control,Project.State);

%% Baseline
Out.Baseline.Elec = zeros(8760,1);
Out.Baseline.Heat = zeros(8760,1);
steps = length(demand.DemandE);
Ts = 8760/steps;
for i = 1:1:8760
    Out.Baseline.Elec(i) = sum(demand.DemandE(1+1/Ts*(i-1):1/Ts*i)*Ts);
    Out.Baseline.Heat(i) = sum(demand.DemandH(1+1/Ts*(i-1):1/Ts*i)*Ts);
end
monthDays = [0 31 59 90 120 151 181 212 243 273 304 334 365];
for i = 1:1:length(monthDays)-1
    Out.Baseline.Fuel(i) = sum(demand.DemandH(24/Ts*monthDays(i)+1:24/Ts*monthDays(i+1)))*Ts;
end

%% Dispatch
Out.Dispatch.Elec = Out.eOut.Grid_purchases_hourly_kWh;
monthDays = [0 31 59 90 120 151 181 212 243 273 304 334 365];
for i = 1:1:length(monthDays)-1
    Out.Dispatch.Fuel(i) = sum(Out.eOut.Fuel(24/Ts*monthDays(i)+1:24/Ts*monthDays(i+1))); %don't mulltiply by Ts, already converted from kW to kWh
end
Out.Dispatch.ElecTotProd = Out.eOut.Total_Electricity_produced_kWh;

Out.Dispatch.SysSize = Out.eOut.SysSize;
Out.Dispatch.ChillerSize(1) = Out.eOut.ChillerSize(1);
Out.Dispatch.ChillerSize(2) = Out.eOut.ChillerSize(2);
Out.Dispatch.TESsize = Out.eOut.TESsize;
Out.Dispatch.BatterySize = Out.eOut.BatterySize;
Out.Dispatch.RenewableSize(1) = Out.eOut.RenewableSize(1);
Out.Dispatch.RenewableSize(2) = Out.eOut.RenewableSize(2);

%% Financial calculations
Out.costOut=FinancialCalcs2(Out.Baseline,Out.Dispatch,Project.Utilities.Grid,Project.Economic,Project.State,user,Project.Utilities.NatGas,scale);

close(MSGra)