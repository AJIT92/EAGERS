%%CreateTestGenerators
global Model_dir MapWaitbarHandle OnOff Plant DateSim scaleTime Last24hour RealTimeData HistProf Si LB UB CurrentState
OnOff = [];

E = menu('Select Option','Create New Generator Set', 'Create Demand profile','Develop Perfomance Map for Generator Set','Test Generator Dispatch');
if E ==1
    %% load and size generators
    load(fullfile(Model_dir,'component library','GeneratorList'))
    list = {};
    for i= 1:1:length(Generators)
        list(end+1) = cellstr(Generators(i).Name);
    end
    [s,v] = listdlg('PromptString','Select Distributed Energy Resources to Consider','SelectionMode','multiple','ListString',list);
    Generators = Generators(s);
    list = {};
    sizes = {};
    for i= 1:1:length(Generators)
        list(end+1) = cellstr(Generators(i).Name);
        sizes(end+1) = cellstr(num2str(Generators(i).Size));
    end
    %select sizes
    s = str2double(inputdlg(list,'Specify Size in kW or kWh',1,sizes));
    for i= 1:1:length(Generators)
        if strcmp(Generators(i).Type,'Electric Generator') || strcmp(Generators(i).Type,'CHP Generator')||strcmp(Generators(i).Type,'Chiller')%electric or CHP gen
            Generators(i).VariableStruct.StartCost = Generators(i).VariableStruct.StartCost*s(i)/Generators(i).Size;
        end
        Generators(i).Size = s(i);
        if strcmp(Generators(i).Type,'Electric Generator') || strcmp(Generators(i).Type,'CHP Generator')%electric or CHP gen
            Generators(i).VariableStruct.Startup.Electricity = linspace(0,Generators(i).Output.Capacity(2),19)*Generators(i).Size;
            Generators(i).VariableStruct.Startup.Input = Generators(i).VariableStruct.Startup.Electricity/Generators(i).Output.Electricity(2);
            Generators(i).VariableStruct.Shutdown.Electricity = linspace(Generators(i).Output.Capacity(2),0,19)*Generators(i).Size;
            Generators(i).VariableStruct.Shutdown.Input = Generators(i).VariableStruct.Shutdown.Electricity/Generators(i).Output.Electricity(2);
            Generators(i).VariableStruct.Startup.Heat = Generators(i).VariableStruct.Startup.Input*Generators(i).Output.Heat(2);
            Generators(i).VariableStruct.Shutdown.Heat = Generators(i).VariableStruct.Shutdown.Input*Generators(i).Output.Heat(2);
        end
        if strcmp(Generators(i).Type,'Chiller') %electric chiller
            Generators(i).VariableStruct.Startup.Cooling = linspace(0,Generators(i).Output.Capacity(2),19)*Generators(i).Size;
            Generators(i).VariableStruct.Startup.Input = Generators(i).VariableStruct.Startup.Cooling./Generators(i).Output.Cooling(2);
            Generators(i).VariableStruct.Shutdown.Cooling = linspace(Generators(i).Output.Capacity(2),0,19)*Generators(i).Size;
            Generators(i).VariableStruct.Shutdown.Input = Generators(i).VariableStruct.Shutdown.Cooling/Generators(i).Output.Cooling(2);
        end
        if strcmp(Generators(i).Type,'Heater') %gas heater
            Generators(i).VariableStruct.Startup.Heat = linspace(0,Generators(i).Output.Capacity(2),19)*Generators(i).Size;
            Generators(i).VariableStruct.Startup.Input = Generators(i).VariableStruct.Startup.Heat./Generators(i).Output.Heat(2);
            Generators(i).VariableStruct.Shutdown.Heat = linspace(Generators(i).Output.Capacity(2),0,19)*Generators(i).Size;
            Generators(i).VariableStruct.Shutdown.Input = Generators(i).VariableStruct.Shutdown.Heat/Generators(i).Output.Heat(2);
        end
    end
    f=char(inputdlg('Name','Specify Name For Generator Set',1,{'genset'}));
    f = fullfile(Model_dir,'Plant','GenSets',f);
    save(f,'Generators')
end
if E ==2 %create a load profile
    %%do this manually, create . mat file with demands normalized to 1.
    
end
if E ==3 || E==4
    files=dir(fullfile(Model_dir,'Plant','GenSets','*.mat'));
    list=strrep({files.name},'.mat','');
    [s,v] = listdlg('PromptString','Select Generator Set','SelectionMode','single','ListString',list);
    load(char(list(s)));
    Plant.Generator = Generators;
end

%% Load Demand
DateSim = 736330;%this is Jan 1 in 2016
Plant.Data.Timestamp(1) = DateSim;
if E==4
    files=dir(fullfile(Model_dir,'Plant','LoadProfiles','*.mat'));
    list=strrep({files.name},'.mat','');
    [s,v] = listdlg('PromptString','Select Load Profile','SelectionMode','single','ListString',list);
    load(char(fullfile('LoadProfiles',list(s))));
    s = fields(Demand);
    list = {'Duration (hours)'; 'Resoultion (.25 = 15min)';};
    defaultAns = {'24';'1'};
    Q = {'H', 'Heating','200'; 'E', 'Electric','400'; 'C', 'Cooling','400';};
    for i = 1:1:length(s)
        a = find(strcmp(Q(:,1),s(i)));
        list(end+1) = cellstr(strcat('Peak Demand: ',Q(a,2)));
        defaultAns(end+1) = cellstr(Q(a,3));
    end
    m = str2double(inputdlg(list,'Specify Load',1,defaultAns));
    for i = 1:1:length(s)
        Demand.(char(s(i))) = m(2+i)*Demand.(char(s(i)))(4*m(2):4*m(2):4*m(1))';
    end
elseif E ==3
    m=100;
    clear Demand
    for i = 1:1:length(Generators)
        if strcmp(Generators(i).Type,'CHP Generator') ||strcmp(Generators(i).Type,'Electric Generator')
            Demand.E = zeros(1,m);
        elseif strcmp(Generators(i).Type,'Chiller')%% Cooling Generators 
            Demand.C = zeros(1,m);
        elseif strcmp(Generators(i).Type,'Heater') || strcmp(Generators(i).Type,'CHP Generator')%% Heaters 
            Demand.H = zeros(1,m);
        end
    end
end
%% Specify Electric utility & gas if necessary
if (E ==3 || E==4) && isfield(Demand,'E')
    n = length(Generators)+1;
    files=dir(fullfile(Model_dir,'component library','Utility','*.mat'));
    list=strrep({files.name},'.mat','');
    [s,v] = listdlg('PromptString','Select Utility Profile','SelectionMode','single','ListString',list);
    load(char(fullfile('Utility',list(s))))
    
    Plant.Generator(n).Type = component.Type;
    Plant.Generator(n).Name = component.Name;
    Plant.Generator(n).Source = 'Electricity';
    Plant.Generator(n).Output = [];
    Plant.Generator(n).Size = inf;
    Plant.Generator(n).Enabled = 1;
    Plant.Generator(n).VariableStruct = component;
    Plant.Generator(n).VariableStruct.MinImportThresh = 0;

    s = str2double(inputdlg({'A gas utility may be needed. Gas Rate in $ per MMBTU';},'Gas Rate',1,{'5.00'}));
    Plant.Generator(n+1).Type = 'Utility';
    Plant.Generator(n+1).Name = 'Gas Supply';
    Plant.Generator(n+1).Source = 'NG';
    Plant.Generator(n+1).Enabled = 1;
    Plant.Generator(n+1).VariableStruct.WinRateTable = ones(size(Plant.Generator(n).VariableStruct.WinRateTable));
    Plant.Generator(n+1).VariableStruct.SumRateTable = ones(size(Plant.Generator(n).VariableStruct.SumRateTable));
    Plant.Generator(n+1).VariableStruct.WinRate = s;
    Plant.Generator(n+1).VariableStruct.SumRate = s;
    Plant.Generator(n+1).VariableStruct.Timestamp = (0:365*24+1)/24+DateSim;
    Plant.Generator(n+1).VariableStruct.Rate = s*ones(365*24+2,1);
end


%% 
if E == 3%% Start iteration, solving quadratic problem at loads in the range minGen to maxGen
    Resolution = 1;
    Plant.optimoptions.Horizon = Resolution;
    Plant.optimoptions.tspacing = 'constant';
    Plant.optimoptions.Resolution = Resolution;
    Plant.optimoptions.Topt = 3600*Resolution;
    Plant.optimoptions.Tmpc = 3600*Resolution;
    Plant.optimoptions.sequential = 0;
    Plant.optimoptions.ExcessHeat = 1;
    Plant.optimoptions.thresholdSteps = 1;
    Time = buildTimeVector(Plant.optimoptions);
    scaleTime = 1;
    loadGenerator 
    %% load demands
    if isempty(Plant.Generator(n).VariableStruct.MinImportThresh)
        Esize = 0;
        EminSize = [];
    else Esize = Plant.Generator(n).VariableStruct.MinImportThresh;
        EminSize = Plant.Generator(n).VariableStruct.MinImportThresh;
    end
    Csize = 0;
    CminSize = [];
    Hsize = 0;
    for i = 1:1:length(Generators)
        if strcmp(Generators(i).Type, 'Electric Generator') ||strcmp(Generators(i).Type,'CHP Generator')
            Esize = Esize + Generators(i).Size;
            EminSize(end+1) = Generators(i).VariableStruct.Startup.Electricity(end);
            if strcmp(Generators(i).Type,'CHP Generator')
                Hsize = Hsize + Generators(i).Size*Generators(i).Output.Heat(end)/Generators(i).Output.Electricity(end);
            end
        elseif strcmp(Generators(i).Type, 'Chiller')
            Csize = Csize + Generators(i).Size;
            CminSize(end+1) = Generators(i).VariableStruct.Startup.Cooling(end);
        elseif strcmp(Generators(i).Type, 'Heater')
            Hsize = Hsize + Generators(i).Size;
        end

    end
    EminSize = min(EminSize);
    CminSize = min(CminSize);
    if isfield(Demand,'E')
        Demand.E = linspace(EminSize,Esize,m);
        if isfield(Demand,'C')
            s = str2double(inputdlg(strcat('Specify the cooling demand (kW): max size is ',num2str(Csize)),'Cooling Demand', 1,{num2str(Csize)}));
            Demand.C = s + Demand.C;
        end
        if isfield(Demand,'H')
            s = str2double(inputdlg(strcat('Specify the heating demand (kW): max size is ',num2str(Hsize)),'Heating Demand', 1,{num2str(Hsize)}));
            Demand.H = s+Demand.H;
        end
    elseif isfield(Demand,'C')
        Demand.C = linspace(CminSize,Csize,m);
    end
    tic
    scaleCost = updateGeneratorCost(0);
    marginCost = updateMarginalCost(ones(nS+1,1)*UB,scaleCost,Time);%assume everything on to find a marginal cost.
    MapWaitbarHandle=waitbar(0,'Building Optimal Dispatch Map');
    [costMin, OptimalOutput] = StepByStepDispatch(Demand,scaleCost,ones(m,1)*Resolution,[],[],marginCost,[]);%an empty limit denotes unconstrained
    close(MapWaitbarHandle)
    MapWaitbarHandle =[];
    plotDispatchMap(Generators,OptimalOutput,costMin,Demand)
    toc
end

if E == 4% Test Dynamic Dispatch
    Horizon= m(1);
    Resolution = m(2);
    %%
    Plant.optimoptions.Horizon = Horizon;
    Plant.optimoptions.tspacing = 'constant';
    Plant.optimoptions.Resolution = Resolution;
    Plant.optimoptions.Topt = 3600;
    Plant.optimoptions.Tmpc = 600;
    Plant.optimoptions.ExcessHeat = 1;
    Plant.optimoptions.thresholdSteps = 6;
    Plant.optimoptionsBuffer = 0; % percentage for buffer on storage
    Time = buildTimeVector(Plant.optimoptions);
    Outs = fieldnames(Demand);
    if ismember('C',Outs) && ismember('E',Outs)
        C = menu('Select Chiller Dispatch Option','Dispatch Chillers 1st, then generators', 'Dispatch Chillers with Generators');
        if C==2
            Plant.optimoptions.sequential = 0;
            Outs = {'E'};
        else
            Plant.optimoptions.sequential = 1;
        end
    else
        Plant.optimoptions.sequential = 1;
    end
    %% Load Historicals
    Si = 2;
    nS = length(Time);
    RealTimeData.Timestamp = [0,Time]./24+DateSim;
    RealTimeData.Temperature = ones(1,length(RealTimeData.Timestamp));
    HistProf.Timestamp = Time./24+DateSim;
    HistProf.Temperature= ones(12,nS);
    %create the surface of the demand so that updateMatrices forecasts the
    %exact demand that is input
    x = zeros(nS*2+2,1);
    y = zeros(nS*2+2,1);
    x(1:2:nS*2+1,:) = [0,Time];
    x(2:2:end) = [0,Time];
    y(1:2:end) = 2;%y is the temperature variable and goes from 0 to 2, temp is always 1 in TestDisp
    z = zeros(nS*2+2,1);
    S = fields(Demand);
    for i = 1:1:length(S)
        ID.(char(S(i))) = Demand.(char(S(i)))(1); %Initial Demand
        dem = [ID.(char(S(i))), Demand.(char(S(i)))]';
        z(1:2:end) = dem;
        z(2:2:end) = dem;
        RealTimeData.Demand.(char(S(i))) = [ID.(char(S(i))), Demand.(char(S(i)))];
        HistProf.(char(S(i))) = fit([x,y],z,'linearinterp');
    end
    nG = length(Plant.Generator);
    
    %% Load generators and find initial condition
    scaleTime = 1;
    loadGenerator 
    stor = zeros(1,nG);
    for i=1:1:nG
        if isfield(Plant.Generator(i).OpMatA,'Stor')
            stor(i) = i;
        end
    end
    scaleCost = updateGeneratorCost(Time);
    Data_t0 = GetCurrentData(DateSim);
    IC = StepByStepDispatch(Data_t0.Demand,scaleCost(1,:),Resolution, [],'',[]);
    IC(stor>0) = .5*UB(stor>0); % IC = halfway charged energy storage
    OnOff = logical(IC>=LB);
    CurrentState.Generators=IC;
    nS_step1 = round(Horizon/Time(1));
    Last24hour =[];%eliminate any old data stored here
    %need to have this in terms of the first timestep
    Last24hour = updateForecast(DateSim-1,(1:nS_step1).*Time(1));%need to have one extra timestep for prediction
    Last24hour.Timestamp = DateSim-1+((0:nS_step1-1).*(Time(1)/24));
    
    %% load fitA, update matrices and run 1st step
    Organize = Plant.OpMatA.Organize;
    marginCost = updateMarginalCost(UB,scaleCost(1,:),Time);
    [QPall,~] = updateMatrices(Plant.OpMatA.QP,Plant.OpMatA.Organize,IC,Time,scaleCost,marginCost,[]); 
    for i = 1:1:length(stor)% add extra constraint for end of storage to be IC of storage Currently does not seperate storage by category if more than 1 QPall (i.e. 'C' and 'E' sequential)
        if stor(i)~=0
            k = Organize.States{stor(i)};
            k = k(nS+1); %end condition
            QPall.(Outs{1}).Aeq(end+1,k) = 1;
            QPall.(Outs{1}).beq(end+1) = IC(stor(i));
        end
    end
    [FirstDisp, ~, Feasible] = DispatchQP(QPall,Organize,[]);%this is the dispatch with fit A
    
    %% use modified generator demand profile for Step 2
    Out = Plant.optimoptions.Outputs;
    %change demand to not include the amount covered by storage
    Demand = AccountForSelfDischarge(Demand,Time);
    for j = 1:1:length(Out);
        GenDemand.(Out{j}) = zeros(nS,1);
        for i = 1:1:length(Plant.Generator)
            if isfield(Plant.Generator(i).OpMatB.output,Out{j})
                if ~isfield(Plant.Generator(i).OpMatB,'Stor') %don't include the storage use/charging in this new profile;
                    GenDemand.(Out{j}) = GenDemand.(Out{j}) + FirstDisp(2:end,i)*Plant.Generator(i).OpMatB.output.(Out{j});
                end
            end
        end
    end
    tic
    OptimalState = StepByStepDispatch(GenDemand,scaleCost,ones(nS,1)*Plant.optimoptions.Resolution,IC, 'initially constrained',FirstDisp);
    disp(['time for step 2 is ' num2str(toc)])
    %% check if any gen can be completley removed
    Locked = true(nS+1,nG);
    %% Rule 1 if off-line for entire optimal dispatch, then the generator should be off
    redoStep1 = false;
    for i = 1:1:nG
        if nnz(OptimalState(:,i))==0
            Locked(:,i)= false;
            redoStep1 = true;
        end
    end
    if redoStep1
        [FirstDisp, ~, Feasible] = DispatchQP(QPall,Organize,Locked);%this is the dispatch with fit A
        if Feasible ~=1
            Locked = true(nS+1,nG); %reset to all on default
        end
    end
    %% load & run fitB for step 3
    tic
    Organize = Plant.OpMatB.Organize;
    marginCost = updateMarginalCost(FirstDisp,scaleCost,Time);
    [QPall,~] = updateMatrices(Plant.OpMatB.QP,Organize,IC,Time,scaleCost,marginCost,[]); %update fit B matrices
    for i = 1:1:length(stor)% add extra constraint for end of storage to be IC of storage Currently does not seperate storage by category if more than 1 QPall (i.e. 'C' and 'E' sequential)
        if stor(i)~=0
            k = Organize.States{stor(i)};
            k = k(nS+1); %end condition
            QPall.(Outs{1}).Aeq(end+1,k) = 1;
            QPall.(Outs{1}).beq(end+1) = IC(stor(i));
        end
    end
    [GenDisp, dX] = FilterGenerators(QPall,Organize,IC,Demand,FirstDisp,OptimalState,scaleCost,Locked);
    disp(['time to load, update matrices and run step 3 is ' num2str(toc)])
    
    %% record input to asses cost of dispatch then plot
    Input = 0*GenDisp;
    chillers = zeros(1, length(nG));
    for i = 1:1:nG
        if ~isempty(Plant.Generator(i).Output)
            cap = Plant.Generator(i).Output.Capacity*UB(i);
        end
        eff = [];
        if strcmp(Plant.Generator(i).Type,'Electric Generator') || strcmp(Plant.Generator(i).Type,'CHP Generator')
            eff = Plant.Generator(i).Output.Electricity;
        elseif strcmp(Plant.Generator(i).Type,'Chiller') &&  ~ismember('E',Outs)%don't include cost if it shows up in generator demand
            eff = Plant.Generator(i).Output.Cooling;
            chillers(i) = i;
        elseif strcmp(Plant.Generator(i).Type,'Heater')
            eff = Plant.Generator(i).Output.Heat;    
        end
        if ~isempty(eff)
            Input(:,i) = GenDisp(:,i)./interp1(cap,eff,GenDisp(:,i));
        elseif strcmp(Plant.Generator(i).Type,'Utility')
            Input(:,i) = GenDisp(:,i);
        end
    end
    Input(isnan(Input))=0;
    if ~Plant.optimoptions.sequential% if you are running chiller and generator at the same time, then the chiller electric use is added to demand, so don't count chiller use twice
        Input(:,chillers>0) = 0;
    end
    Cost = sum(NetCostCalc(GenDisp,Time,Input));
    disp(['Optimal dispatch cost is : ', num2str(Cost)])
    plotDispatch2(GenDisp,Demand,Time,1)
    
    %% Run Baseline% will keep all generators on at all times. May not be feasible
    MapWaitbarHandle=waitbar(0,'Running Baseline Dynamic Economic Dispatch');
    Organize = Plant.OpMatB.Organize;
    marginCost = updateMarginalCost(UB,scaleCost(1,:),Time);
    ICallOn = LB;
    ICallOn(stor>0) = .5*UB(stor>0); % IC = halfway charged energy storage
    OnOff = logical(IC>=LB);
    CurrentState.Generators=IC;
    Locked = [ICallOn>=LB;ones(length(Time),length(UB))]; 
    [QPall,~] = updateMatrices(Plant.OpMatB.QP,Organize,ICallOn,Time,scaleCost,marginCost,[]); %update fit B matrices starting with all gen @ LB
    for i = 1:1:length(stor)% add extra constraint for end of storage to be IC of storage Currently does not seperate storage by category if more than 1 QPall (i.e. 'C' and 'E' sequential)
        if stor(i)~=0
            k = Organize.States{stor(i)};
            k = k(nS+1); %end condition
            QPall.(Outs{1}).Aeq(end+1,k) = 1;
            QPall.(Outs{1}).beq(end+1) = IC(stor(i));
        end
    end
    [Baseline.GeneratorState,  ~, Feasible] = DispatchQP(QPall,Organize,Locked);%this is the dispatch with fit B, and all generators on
    close(MapWaitbarHandle)
    MapWaitbarHandle =[];
    %% asses costs and plot baseline
    if Feasible==1
        Input = 0*Baseline.GeneratorState;
        for i = 1:1:nG
            if ~isempty(Plant.Generator(i).Output)
                cap = Plant.Generator(i).Output.Capacity*UB(i);
            end
            eff = [];
            if strcmp(Plant.Generator(i).Type,'Electric Generator') || strcmp(Plant.Generator(i).Type,'CHP Generator')
                eff = Plant.Generator(i).Output.Electricity;
            elseif strcmp(Plant.Generator(i).Type,'Chiller') && ~ismember('E',Outs)%don't include cost if it shows up in generator demand
                eff = Plant.Generator(i).Output.Cooling;
            elseif strcmp(Plant.Generator(i).Type,'Heater')
                eff = Plant.Generator(i).Output.Heat;    
            end
            if ~isempty(eff)
                Input(:,i) = Baseline.GeneratorState(:,i)./interp1(cap,eff,Baseline.GeneratorState(:,i));
            elseif strcmp(Plant.Generator(i).Type,'Utility')
                Input(:,i) = GenDisp(:,i);
            end
        end
        Input(isnan(Input))=0;
        Baseline.NetCost = NetCostCalc(Baseline.GeneratorState,Time,Input);
        disp(['Standard Dynamic Economic Dispatch cost is : ', num2str(sum(Baseline.NetCost))])
        plotDispatch2(Baseline.GeneratorState,Demand,Time,10)
    else disp('Baseline simulation is infeasible');
    end
end