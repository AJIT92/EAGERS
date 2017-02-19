function dX = RunLinSystem(t,X)
%Finds dX's of a system with LinMod's array of LTI models
global Tags LinMod IterCount TagInf WaitBar SimSettings Jcount
if isempty(TagInf) %create TagInf Fields
    TagInf.Time(1) =0;
    Jcount = 0;
end
if t==0 || t == TagInf.Time(IterCount)
    Jcount = Jcount+1;
elseif t>TagInf.Time(IterCount)%+1e-5
    IterCount =IterCount+1;
    Jcount = length(X);
elseif IterCount>1 && t<TagInf.Time(IterCount)
    Jcount = 1;
    while t<TagInf.Time(IterCount-1)
        IterCount = IterCount-1;
    end
end
TagInf.Time(IterCount,1) = t;

InterpVal = X(1); 
c= [0 0 0 0];
n = length(LinMod.InterpVec);
if LinMod.InterpVec(1) > LinMod.InterpVec(2) %values of InterpVec are decreasing
    k = nnz(LinMod.InterpVec > InterpVal);
else
    k = nnz(LinMod.InterpVec < InterpVal);
end
if k == n
    nModel(4) = n;
    c(4) = 1; 
elseif k==0
    nModel(1) = 1;
    c(1) = 1; 
elseif k == n-1
    nModel(4) = k+1;
    nModel(3) = k;
    c(4) = (InterpVal - LinMod.InterpVec(n-1))/(LinMod.InterpVec(n) - LinMod.InterpVec(n-1));
    c(3) = (InterpVal - LinMod.InterpVec(n))/(LinMod.InterpVec(n-1) - LinMod.InterpVec(n));
elseif k == 1
    nModel(2) = k+1;
    nModel(1) = k;
    c(2) = (InterpVal - LinMod.InterpVec(1))/(LinMod.InterpVec(2) - LinMod.InterpVec(1));
    c(1) = (InterpVal - LinMod.InterpVec(2))/(LinMod.InterpVec(1) - LinMod.InterpVec(2));
else 
    nModel(4) = k+2;
    nModel(3) = k+1;
    nModel(2) = k;
    nModel(1) = k-1;
    Vec = [LinMod.InterpVec(k-1),LinMod.InterpVec(k),LinMod.InterpVec(k+1),LinMod.InterpVec(k+2)];
    c = ones(1,4);
    for i = 1:4
        for j = 1:4
            if i ~= j
                c(i) = c(i)*(InterpVal - Vec(j))/(Vec(i) - Vec(j));
            end
        end
    end
end

%% iterate linear model + controller to find dX and dUX: I would like to find a way for this iteration to converge faster
nS = length(LinMod.Model{1}.A);%number of model states
ControlStates = X(nS+1:end,1); %seperate the controller states (at end of vector)
X = X(1:nS,1); %reduce to just the model states

U = Tags.U; %controller outputs
O = zeros(length(LinMod.Model{1}.Out0),1); %model Out
dX = zeros(length(X),1); %change in plant states
dUX = zeros(length(ControlStates),1); %change in controller states

repeat = true;
count = 0;
while repeat
    U_old = U;
    O_old = O;
    dX_old = dX;
    dUX_old = dUX;
    dX = 0;
    O = 0;
    for i = 1:1:4
        if c(i)>0
            dX = dX + c(i)*(LinMod.Model{nModel(i)}.A*(X-LinMod.Model{nModel(i)}.X0) + LinMod.Model{nModel(i)}.B*(U-LinMod.Model{nModel(i)}.U0));%calculate dY from Y's and U's  %% Issue coming from Y causing dY's
            O = O + c(i)*(LinMod.Model{nModel(i)}.Out0 + (LinMod.Model{nModel(i)}.C*(X-LinMod.Model{nModel(i)}.X0) + LinMod.Model{nModel(i)}.D*(U-LinMod.Model{nModel(i)}.U0)));%calculate O from Y's and U's
        end
    end
    [U,dUX] = Controller(t,O,ControlStates);
    relTol = abs([O;U;dX;dUX;]-[O_old;U_old;dX_old;dUX_old;])./abs([O;U;dX;dUX;]);
    relTol(abs([O;U;dX;dUX;])<1e-14) = [];
    if ~any(relTol>1e-3) %relative tolerance
        repeat = false;
    end
    count = count+1;
end
% TagInf.U(IterCount,:) = U;
% TagInf.O(IterCount,:) = O;
Tags.U = U;
dX(end+1:end+length(ControlStates)) = dUX;

if t>0 && WaitBar.Show == 1 && Jcount==length(X) && isfield(LinMod,'Scope')
    n = length(LinMod.Scope);
    dt = TagInf.Time(IterCount,1) - TagInf.Time(IterCount-1,1);
    points = nnz(TagInf.Time>(TagInf.Time(IterCount,1)-100*dt));
    for i = 1:1:n
        tagName = LinMod.Scope{i};
        r = strfind(tagName,'.');
        block = tagName(1:r-1);
        tag = tagName(r+1:end);
        if isfield(TagInf,block) && isfield(TagInf.(block),tag)
            figure(i)
            plot(TagInf.Time(IterCount-points+1:IterCount,1),TagInf.(block).(tag)(IterCount-points+1:IterCount,:))
            ylabel(tagName);
            xlabel('Time in seconds');
        end
    end
end
%% update waitbar
if t ==0
    Text = 'calculating Jacobian';
    x = Jcount/length(X);
elseif IterCount>1 && Jcount<length(X)
    Text = 're-calculating Jacobian';
    x = Jcount/length(X);
else
    if strcmp(WaitBar.Text,'Running linear model with transient')
        x = t/SimSettings.RunTime;
    else x = max(0,(log(t/1e-4))/(log(SimSettings.RunTime/1e-4)));
    end
    Text = {WaitBar.Text;strcat('    Time = ', num2str(t))};
end
if WaitBar.Show == 1
    waitbar(x,WaitBar.Handle,Text);
end

function [U,dUx] = Controller(t,ModelOut,ControlStates)
% what to do if the controller outputs and # of states don't align?
global LinMod Tags TagInf IterCount

n=0; %counter for inlet values
m = 0; %counter for outlet values
list = fieldnames(Tags.ModelInput);
for b = 1:1:length(list)
    block = list{b};
    ports = fieldnames(Tags.ModelOutput.(block)); %outputs of model (inputs to controller)
    for j = 1:1:length(ports)
        port = ports{j};
        if isstruct(Tags.ModelOutput.(block).(port))
            f = fieldnames(Tags.ModelOutput.(block).(port));
            for k = 1:1:length(f)
                s = Tags.ModelOutput.(block).(port).(f{k});
                Inlet.(port).(f{k}) = ModelOut(n+1:n+length(s));
                n = n+length(s);
            end
        else
            s = Tags.ModelOutput.(block).(port);
            Inlet.(port) = ModelOut(n+1:n+length(s));
            n = n+length(s);
        end
    end
    Out = feval(LinMod.Controls.(block).type,t,ControlStates,Inlet,LinMod.Controls.(block),'Outlet');
    ports = fieldnames(Tags.ModelInput.(block)); %inputs to model (outputs from controller)
    for j = 1:1:length(ports)
        port = ports{j};
        if isstruct(Out.(port))
            f = fieldnames(Out.(port));
            for k = 1:1:length(f)
                s = Out.(port).(f{k});
                U(m+1:m+length(s),1) = s;
                m = m+length(s);
            end
        else
            s = Out.(port);
            U(m+1:m+length(s),1) = s;
            m = m+length(s);
        end
    end
    dUx = feval(LinMod.Controls.(block).type,t,ControlStates,Inlet,LinMod.Controls.(block),'dY');
    
    if isfield(LinMod.Controls.(block),'TagInf')
        tagNames = LinMod.Controls.(block).TagInf;
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
end