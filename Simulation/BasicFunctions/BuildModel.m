function BuildModel(Plant)
%builds a single model from all of the component blocks, makes connections
%between inlets and outlets
global modelParam SimSettings Inlet Outlet TagInf TagFinal IterCount WaitBar Tags
modelParam=[]; Inlet=[]; Outlet =[]; 
tic; 
if ~isfield(WaitBar,'Show')
    WaitBar.Show = 1;
end

if ~isfield(Tags,'Options')
    Tags.Options.ODESet = [];
end

n = 0; %counter for the states
if isfield(Plant,'Scope')
    modelParam.Scope = Plant.Scope;
end
if isfield(Plant,'Plot')
    modelParam.Plot = Plant.Plot;
end
if isfield(Plant,'NominalPower')
    modelParam.NominalPower = Plant.NominalPower;
end
modelParam.Components = Plant.Components;
modelParam.Controls = Plant.Controls;
CompNames = fieldnames(modelParam.Components);
controls = fieldnames(modelParam.Controls);
list = [CompNames;controls;];

%% Load block initializations
for k = 1:1:length(list)
    block = list{k};
    if any(strcmp(controls,block))
        Co = 'Controls';
    else Co = 'Components';
    end
    F = strcat('Initialize',modelParam.(Co).(block).type);
    modelParam.(block) = feval(F,modelParam.(Co).(block));
    states = modelParam.(block).IC; %initial condition
    scale = modelParam.(block).Scale; %scalar on states so they are all ~ 1
    s = length(states);
    modelParam.(block).States = n+1:n+s;
    if s>0
        modelParam.IC(n+1:n+s,1) = states; 
        modelParam.Scale(n+1:n+s,1) = scale; 
        n = n+s;
    end 
end
connectPorts
%% Converge component initializations to an approximation of steady-state operation
SimSettings.RunTime = 3600*24;
if isfield(modelParam,'NominalPower')
    SimSettings.PowerDemand    = [1 1]*modelParam.NominalPower;
    SimSettings.PowerTime = linspace(0,SimSettings.RunTime,length(SimSettings.PowerDemand));
end
OldInlet = Inlet; %initial condition inlets
Error = 1;
error = zeros(length(controls),1);
count = 0;
Tol = 1e-2;
WaitBar.Text = {'Converging Initial Guess';'    ratio error:tolerance = inf'};
WaitBar.Handle=waitbar(0,WaitBar.Text);
while abs(Error)>Tol %repeat several times to propogate inlets to outlets
    SolvePinitial(list)%% Solve for initial Pressures 
    for k = 1:1:length(list)
        block = list{k};
        if any(strcmp(controls,block))
            Co = 'Controls';
        else Co = 'Components';
        end           
        Inlet.(block) = RefreshInlet(block);
        if HasInletChanged(Inlet.(block),OldInlet.(block))
            F = strcat('Initialize',modelParam.(Co).(block).type);
            modelParam.(block) = feval(F,modelParam.(block),Inlet.(block));
            newStates = length(modelParam.(block).Scale)-length(modelParam.(block).States);
            if newStates ==0
                modelParam.Scale(modelParam.(block).States,1) = modelParam.(block).Scale;
                modelParam.IC(modelParam.(block).States,1) = modelParam.(block).IC;
            else
                PrevStates = modelParam.Scale(1:modelParam.(block).States(1)-1,1);
                AfterStates = modelParam.Scale(modelParam.(block).States(end)+1:end,1);
                modelParam.Scale = [PrevStates; modelParam.(block).Scale;AfterStates;];
                PrevIC = modelParam.IC(1:modelParam.(block).States(1)-1,1);
                AfterIC = modelParam.IC(modelParam.(block).States(end)+1:end,1);
                modelParam.IC = [PrevIC; modelParam.(block).IC;AfterIC;];
                modelParam.(block).States = modelParam.(block).States(1):(modelParam.(block).States(1)+length(modelParam.(block).Scale)-1);
                for n = k+1:1:length(list) 
                    block2 = list{n};
                    modelParam.(block2).States = modelParam.(block2).States + newStates;
                end
            end
            list2 = modelParam.(block).OutletPorts;
            for i = 1:1:length(list2) %update inlets connected to outlets
                Outlet.(block).(list2{i}) = modelParam.(block).(list2{i}).IC;
            end
            OldInlet.(block) = Inlet.(block); %just ran this block, if any changes are made to inlet, need to re-run
        end
        if ismember(block,controls)
            error(k) = modelParam.(block).InitializeError;
        end
    end
    count = count+1;
    Error = max(abs(error));
    if Error<Tol && count<2
        Error = max(Error,1/2^(2*count+1));
    end
    convergeInlets(modelParam.IC)%ensures all of the flow rates & species propogate through blocks that don't have states for species.
    WaitBar.Text = {'Converging Initial Guess';strcat('    ratio error:tolerance = ', num2str(Error/Tol));};
    x = min(1,max(0,1-log(Error/Tol)/5));
    waitbar(x,WaitBar.Handle,WaitBar.Text);
end
close(WaitBar.Handle); 
n = 0;
%% Plot the intial guess of temperature distribution for fuel cells andelectrolyzers
for i = 1:1:length(CompNames)
    if strcmp(modelParam.Components.(CompNames{i}).type,'FuelCell')|| strcmp(modelParam.Components.(CompNames{i}).type,'Electrolyzer')
        block = modelParam.(CompNames{i});
        CellMap(modelParam.IC',block,n+1);
        title(strcat('Inital temperaturedistribution for :',CompNames{i}))
        n = n+1;
    end
end 

%% Run non-linear model to steady state (24 hours)
options = odeset('RelTol',1e-3);
IterCount = 1; TagInf =[]; TagFinal =[]; WaitBar.Text = 'Initialization'; WaitBar.Handle=waitbar(0,WaitBar.Text);
[T, Y] = ode15s(@RunBlocks, [0, SimSettings.RunTime], modelParam.IC,options); disp(strcat('Initialization:',num2str(toc),' seconds'));close(WaitBar.Handle); 
%     PlotSimulation(T,Y,1,0,1)
modelParam.IC = Y(end,:)';


function Change = HasInletChanged(New,Old)
Change = 0;
list2 = fieldnames(New);
for i = 1:1:length(list2)
    port = list2{i};
    if Change ==0 %skip if we already know it has changed
        if isnumeric(New.(port))
            Change = max(New.(port)~=Old.(port));
        elseif isstruct(New.(port))
            f = fieldnames(New.(port));
            for j = 1:1:length(f)
                if ~isfield(Old.(port),f{j})
                    Change = 1;
                elseif Change ~= 1
                    Change = max(New.(port).(f{j})~=Old.(port).(f{j}));
                end
            end
        end
    end
end

function InletBlock = RefreshInlet(block)
global modelParam Outlet Tags
list = modelParam.(block).InletPorts;
for i = 1:1:length(list)
    port = list{i};
    if ~isempty(modelParam.(block).(port).connected)%inlet connected to another outlet
        BlockPort = modelParam.(block).(port).connected{1};
        r = strfind(BlockPort,'.');
        if ~isempty(r)
            connectedBlock = BlockPort(1:r(1)-1);
            if strcmp(connectedBlock,'Tags')
                connectedBlock = BlockPort(r(1)+1:r(2)-1);
                connectedPort = BlockPort(r(2)+1:end);
                InletBlock.(port) = Tags.(connectedBlock).(connectedPort);
            else
                connectedPort = BlockPort(r+1:end);
                InletBlock.(port) = Outlet.(connectedBlock).(connectedPort);
            end
        else
            InletBlock.(port) = feval(BlockPort,0);%finds value of look-up function at time = 0
        end
    else
        InletBlock.(port) = modelParam.(block).(port).IC; %not connected to an outlet or tag
    end
end

function convergeInlets(IC)
global modelParam Inlet Outlet
modelParam.IC = IC;
%% run through actual components to converge Inlets & Outlets
CompNames = fieldnames(modelParam.Components);
list = CompNames;
OldInlet = Inlet; %initial condition inlets
%%run blocks to find outlet conditions & connect inlets
nComp = length(list);
blockSteady = true(nComp,1);
t = 0;
Co = 'Components';
while any(blockSteady)
    for blockCount = 1:1:nComp 
        block = list{blockCount};
        Inlet.(block) = RefreshInlet(block);
        if HasInletChanged(Inlet.(block),OldInlet.(block))
            OldInlet.(block) = Inlet.(block); %next time only run if the inlets have changed from now.
            Outlet.(block) = feval(modelParam.(Co).(block).type,t,IC(modelParam.(block).States).*modelParam.Scale(modelParam.(block).States),Inlet.(block),modelParam.(block),'Outlet');
            blockSteady(blockCount) = true;
        else
            blockSteady(blockCount) = false;
        end
    end
end