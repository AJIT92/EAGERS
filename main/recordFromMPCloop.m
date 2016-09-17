function recordFromMPCloop(varargin)
% This function acts as a) a historian recording all data vales into the datalog, and
%  b)updates the lower frequency records (Dispatch, and Last24hour) used in the dispatch optimization
global Plant Operation Dispatch CurrentState Last24hour DateSim chargeEff dischEff Si %nMPC  Tmpc  
nG = length(Plant.Generator);
Tmpc = Plant.optimoptions.Tmpc;
Time = buildTimeVector(Plant.optimoptions);
dt1 = Time(1);
nMPC = floor(3600*dt1/Plant.optimoptions.Tmpc);
rIndex = (Si-2)*nMPC+2:(Si-1)*nMPC+1; %cumilative index of run points to save in datalog
[Rgen, renew] = RenewableOutput(DateSim,Operation.Timestamp-Operation.Timestamp(1),'Actual');
if ~isempty(renew)
    Dispatch.Dispatch.GeneratorState(Si,renew) = sum(Rgen,1);
    Operation.GeneratorState(:,renew) = Rgen;
    Dispatch.RunData.GeneratorState(rIndex,renew) = Rgen;%add higher frequency renewable Data
end
%%----%%%%% save to .mat each time
%     name = strcat('Sim',num2str(round(now - datenum(2015,1,1))),'Step',num2str(Si));
%     if RealTime
%         f = fullfile(Model_dir,'results','RealTime',strcat(name,'.mat'));
%     else f = fullfile(Model_dir,'results','Virtualization',strcat(name,'.mat'));
%     end
%     save(f,'Operation')%%save sub step response to .mat file    

Dispatch.RunData.Timestamp(rIndex) = Operation.Timestamp;
Dispatch.Dispatch.Timestamp(Si) = DateSim;
Last24hour.Timestamp = [Last24hour.Timestamp(2:end) Operation.Timestamp(end)];
S = fieldnames(Last24hour);
for i = 1:1:length(S)-1 %dont do this loop for fieldname "Timestamp"
    if strcmp(S(i),'T')
        Dispatch.Dispatch.Temperature(Si) = mean(Operation.Temperature);
        Dispatch.RunData.Temperature(rIndex) = Operation.Temperature;
        Last24hour.T(1:end-1) = Last24hour.T(2:end);
        Last24hour.T(end) = Operation.Temperature(end);
    else
        Dispatch.Dispatch.Demand.(S{i})(Si) = mean(Operation.Demand.(S{i}));
        Dispatch.RunData.Demand.(S{i})(rIndex) = Operation.Demand.(S{i});
        Last24hour.(S{i})(1:end-1) = Last24hour.(S{i})(2:end);
        Last24hour.(S{i})(end) = Dispatch.Dispatch.Demand.(S{i})(Si);
    end
end

for i = 1:1:nG 
    if isfield(Plant.Generator(i).OpMatB,'Stor')
        Charge = Operation.GeneratorState(:,i)<0;
        Discharge = Operation.GeneratorState(:,i)>=0;
        for k = 1:1:length(Charge) %convert operation (power out/in) to SOC
            w = (Si-2)*nMPC+k;
            dS = -Operation.GeneratorState(k,i)*Charge(k)*chargeEff(i)*Tmpc/3600 - Operation.GeneratorState(k,i)*Discharge(k)/dischEff(i)*Tmpc/3600;%net charging in kWh
            Dispatch.RunData.GeneratorState(w+1,i) = Dispatch.RunData.GeneratorState(w,i)+dS;
        end
        Dispatch.Dispatch.GeneratorState(Si,i) = CurrentState.Generators(i); %previous energy stored +- energy dispatched
        Dispatch.Dispatch.GeneratorInput(Si,i) = sum(Operation.GeneratorState(Charge,i)*(Tmpc/3600)); % input or charging of energy storage
    else %generators & utilities
        Dispatch.Dispatch.GeneratorState(Si,i) = sum(Operation.GeneratorState(1:nMPC,i))/nMPC; %average power
        Dispatch.RunData.GeneratorState(rIndex,i) = Operation.GeneratorState(:,i);
        Dispatch.Dispatch.GeneratorInput(Si,i) = sum(Operation.GeneratorInput(1:nMPC,i))/nMPC; %average fuel
    end
end