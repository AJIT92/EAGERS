function IC = automaticInitialCondition(Data_t0)
global Plant DateSim
nG = length(Plant.Generator);
scaleCost = updateGeneratorCost(DateSim);%% All costs were assumed to be 1 when building matrices, update Generator costs for the given time
Data_t0.Demand.Renewable = zeros(1,length(Plant.Generator));
for i = 1:1:length(Plant.Generator)
    if strcmp(Plant.Generator(i).Source,'Renewable') && Plant.Generator(i).Enabled
        Data_t0.Demand.Renewable(i) = RenewableOutput(Plant.Generator(i).VariableStruct,DateSim,0,'Actual');
    end
end
IC = StepByStepDispatch(Data_t0.Demand,scaleCost,Plant.optimoptions.Resolution,[],'',[]);
for i=1:1:nG
    if isfield(Plant.Generator(i).OpMatA,'Stor') && Plant.Generator(i).Enabled
        IC(i) = 0.5*Plant.Generator(i).OpMatA.Stor.UsableSize; % IC = halfway charged energy storage
    elseif isfield(Plant.Generator(i).OpMatA,'Stor')
        IC(i) = 0; %storage that is disabled has 0 IC
    end
end

%% specify initial river flow and spillway flows: IC (nL)
if isfield(Plant.Network,'Hydro')
    NodeNames = cell(length(Plant.Network),1);
    for i = 1:1:length(Plant.Network)
        NodeNames(i) = Plant.Network(i).name;
    end
    networkNames = fieldnames(Plant.Network);
    networkNames = networkNames(~strcmp('name',networkNames));
    networkNames = networkNames(~strcmp('Equipment',networkNames));
    for net = 1:1:length(networkNames)
        nLinet(net) = length(Plant.subNet.lineNames.(networkNames{net}));
    end
    nLcum = 0; %cumulative line #
    for net = 1:1:length(networkNames)
        if ~strcmp(networkNames{net},'Hydro')
            nLcum = nLcum+nLinet(net); %all hydro lines have initial state
        else
            for i = 1:1:nLinet(net) 
                IC(nG+nLcum+i) = getHydroFlows(DateSim,i);
            end
        end
    end
end