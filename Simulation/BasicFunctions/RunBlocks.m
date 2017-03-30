function dY = RunBlocks(t,Y)
global modelParam IterCount Inlet Outlet TagInf TagFinal Tags SimSettings WaitBar Jcount

Y = Y.*modelParam.Scale;
if ~isfield(Tags,'Options');
    Tags.Options = [];
end
CompNames = fieldnames(modelParam.Components);
controls = fieldnames(modelParam.Controls);
list = [CompNames;controls;];

if isempty(TagInf) %create TagInf Fields
    TagInf.Time(1) =0;
    Jcount = 0;
    for k = 1:1:length(list)
        block = list{k};
        if any(strcmp(controls,block))
            Co = 'Controls';
        else Co = 'Components';
        end
        if isfield(modelParam.(Co).(block),'TagInf')
            f = modelParam.(Co).(block).TagInf;
            for i = 1:1:length(f)
                TagInf.(block).(f{i}) = [];
            end
        end
    end
end
if t==0 || t == TagInf.Time(IterCount)
    Jcount = Jcount+1;
elseif t>TagInf.Time(IterCount)%+1e-5
    IterCount =IterCount+1;
    Jcount = length(Y);
elseif IterCount>1 && t<TagInf.Time(IterCount)
    Jcount = 1;
    while t<TagInf.Time(IterCount-1)
        IterCount = IterCount-1;
    end
end
TagInf.Time(IterCount,1) = t;

OldInlet = [];%Inlet; %initial condition inlets
%% Update pressure states & inlets
for i = 1:1:length(modelParam.Pstates)
    Pnew = Y(modelParam.Pstates(i));
    if ~isempty(modelParam.Poutlets{i})
        block = modelParam.Poutlets{i}{1};
        port = modelParam.Poutlets{i}{2};
        Outlet.(block).(port) = Pnew; %update outlets based on calculated pressure
    end
end
%% run blocks to find outlet conditions & connect inlets
nComp = length(list);
blockSteady = true(nComp,1);
while any(blockSteady)
    for blockCount = 1:1:nComp 
        block = list{blockCount};
        Inlet.(block) = RefreshInlet(block);
        if ~isfield(OldInlet,block)
            OldInlet.(block) = [];
        end
        if HasInletChanged(Inlet.(block),OldInlet.(block))
            OldInlet.(block) = Inlet.(block); %next time only run if the inlets have changed from now.
            if any(strcmp(controls,block))
                Co = 'Controls';
            else Co = 'Components';
            end
            Outlet.(block) = feval(modelParam.(Co).(block).type,t,Y(modelParam.(block).States),Inlet.(block),modelParam.(block),'Outlet');
            blockSteady(blockCount) = true;
        else
            blockSteady(blockCount) = false;
        end
    end
end

%% run blocks to find dY
dY = 0*Y;
list = [CompNames;controls;];
for k = 1:1:length(list) %run components with states %% record All tags
    block = list{k};
    if nnz(strcmp(CompNames,block))>0
        Co = 'Components';
    else Co = 'Controls';
    end
    if ~isempty(modelParam.(block).IC)%blocks with states
        dY(modelParam.(block).States) = feval(modelParam.(Co).(block).type,t,Y(modelParam.(block).States),Inlet.(block),modelParam.(block),'dY');
    end
    if isfield(modelParam.(block),'TagInf')
        tagNames = modelParam.(block).TagInf;
        for i = 1:1:length(tagNames)
            if isnumeric(Tags.(block).(tagNames{i}))
                TagInf.(block).(tagNames{i})(IterCount,:)=Tags.(block).(tagNames{i});
            else
                f = fieldnames(Tags.(block).(tagNames{i}));
                for j = 1:1:length(f)
                    TagInf.(block).(tagNames{i}).(f{j})(IterCount,:)=Tags.(block).(tagNames{i}).(f{j}); 
                end
            end
        end
    end
    if t== SimSettings.RunTime
        if isfield(modelParam.(block),'TagFinal')
            tagNames = modelParam.(block).TagFinal;
            for i = 1:1:length(tagNames)
                if isnumeric(Tags.(block).(tagNames{i}))
                    TagFinal.(block).(tagNames{i})(IterCount,:)=Tags.(block).(tagNames{i});
                else
                    f = fieldnames(Tags.(block).(tagNames{i}));
                    for j = 1:1:length(f)
                        TagFinal.(block).(tagNames{i}).(f{j})(IterCount,:)=Tags.(block).(tagNames{i}).(f{j}); 
                    end
                end
            end
        end
    end
end
dY = dY./modelParam.Scale;

if t>0 && WaitBar.Show == 1 && Jcount==length(Y) && isfield(modelParam,'Scope')
    n = length(modelParam.Scope);
    dt = TagInf.Time(IterCount,1) - TagInf.Time(IterCount-1,1);
    points = nnz(TagInf.Time>(TagInf.Time(IterCount,1)-100*dt));
    for i = 1:1:n
        tagName = modelParam.Scope{i};
        r = strfind(tagName,'.');
        block = tagName(1:r-1);
        tag = tagName(r+1:end);
        figure(i)
        plot(TagInf.Time(IterCount-points+1:IterCount,1),TagInf.(block).(tag)(IterCount-points+1:IterCount,:))
        ylabel(tagName);
        xlabel('Time in seconds');
    end
end
%% update waitbar
if t ==0
    Text = 'calculating Jacobian';
    x = Jcount/length(Y);
elseif IterCount>1 && Jcount<length(Y)
    Text = 're-calculating Jacobian';
    x = Jcount/length(Y);
else
    if strcmp(WaitBar.Text,'Running non-linear model with transient')
        x = t/SimSettings.RunTime;
    else x = max(0,(log(t/1e-4))/(log(SimSettings.RunTime/1e-4)));
    end
    Text = {WaitBar.Text;strcat('    Time = ', num2str(t))};
end
if WaitBar.Show == 1
    waitbar(x,WaitBar.Handle,Text);
end

function Change = HasInletChanged(New,Old)
Change = false;

if isempty(Old)
    Change = true;
else
    list2 = fieldnames(New);
    for i = 1:1:length(list2)
        port = list2{i};
        if ~Change  %skip if we already know it has changed
            if isnumeric(New.(port))
                Change = comparePort(Old.(port),New.(port));
            elseif isstruct(New.(port))
                f = fieldnames(New.(port));
                for j = 1:1:length(f)
                    if ~isfield(Old.(port),f{j})
                        Change = true;
                    elseif ~Change %skip if we already know it has changed
                        if isnumeric(New.(port).(f{j}))
                            Change = comparePort(Old.(port).(f{j}),New.(port).(f{j}));
                        else
                            g = fieldnames(New.(port).(f{j}));
                            for k = 1:1:length(g)
                                if ~isfield(Old.(port).(f{j}),g{k})
                                    Change = true;
                                elseif ~Change %skip if we already know it has changed
                                    Change = comparePort(Old.(port).(f{j}).(g{k}),New.(port).(f{j}).(g{k}));
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

function C = comparePort(Old,New)
tolerance = eps;%1e-14;
if Old ~= 0
    denom = Old;
else denom = 1;
end
Error = (New - Old)/denom;
if  abs(Error) >= tolerance;
    C = true;
else C = false;
end

function InletBlock = RefreshInlet(block)
global modelParam Outlet Tags
list = modelParam.(block).PortNames;
controls = fieldnames(modelParam.Controls);
for i = 1:1:length(list)
    port = list{i};
    if strcmp(modelParam.(block).(port).type,'in') 
        if ~isempty(modelParam.(block).(port).connected)%inlet connected to another outlet
            BlockPort = char(modelParam.(block).(port).connected);
            r = strfind(BlockPort,'.');
            if ~isempty(r)
                connectedBlock = BlockPort(1:r-1);
                if strcmp(connectedBlock,'Tags')
                    connectedBlock = BlockPort(r(1)+1:r(2)-1);
                    connectedPort = BlockPort(r(2)+1:end);
                    InletBlock.(port) = Tags.(connectedBlock).(connectedPort);
                else
                    connectedPort = BlockPort(r+1:end);
                    if isfield(Outlet.(connectedBlock),connectedPort)
                        InletBlock.(port) = Outlet.(connectedBlock).(connectedPort);
                    else InletBlock.(port) = modelParam.(connectedBlock).(connectedPort).IC;
                    end
                end
            else
                InletBlock.(port) = feval(BlockPort,0);%finds value of look-up function at time = 0
            end
            %%forced interupt between controller and component ports when linearizing model
            if ismember(connectedBlock,controls) && ~ismember(block,controls)
                if isfield(Tags.Options,'AssignedInputs') && Tags.Options.AssignedInputs == 1
                    InletBlock.(port) = Tags.ModelInput.(connectedBlock).(connectedPort);%assign inputs
                else
                    Tags.ModelInput.(connectedBlock).(connectedPort) = InletBlock.(port);%Collect inputs to model (outputs from controller)
                end
            end
            if ismember(block,controls) %&& ~ismember(connectedBlock,fieldnames(modelParam.Controls))
                Tags.ModelOutput.(block).(port) = InletBlock.(port);%collect outputs of model (inputs to controller)
            end
        else
            InletBlock.(port) = modelParam.(block).(port).IC;
        end
    end
end