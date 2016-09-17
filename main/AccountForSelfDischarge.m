function Demand = AccountForSelfDischarge(Demand,Time)
global Plant selfDisch UB dischEff chargeEff
S = fieldnames(Demand);
dt = Time-[0,Time(1:end-1)];
%% account for self-discharge
for i=1:1:length(UB) %Look for any storage types, to add self discharge to demand
    if ismember('E',S) && strcmp(Plant.Generator(i).Type,'Electric Storage')%if battery storage, rewrite to include (self discharge*size)/(discharging eff * charging eff)
        Demand.E = Demand.E + ((selfDisch(i)*UB(i))/(dischEff(i)*chargeEff(i)))*dt;
    elseif ismember('H',S) && strcmp(Plant.Generator(i).Type,'Thermal Storage') && isfield(Plant.Generator(i).OpMatA.output,'H')%if heat storage, rewrite to include (self discharge*size)/(discharging eff * charging eff)
        Demand.H = Demand.H + ((selfDisch(i)*UB(i))/(dischEff(i)*chargeEff(i)))*dt;
    elseif ismember('C',S) && strcmp(Plant.Generator(i).Type,'Thermal Storage') && isfield(Plant.Generator(i).OpMatA.output,'C')%if cold storage, rewrite to include (self discharge*size)/(discharging eff * charging eff)
        Demand.C = Demand.C +  ((selfDisch(i)*UB(i))/(dischEff(i)*chargeEff(i)))*dt;
    end
end