function LinMod = CreateLinModel(SetPoints,HSVtol)
%builds linear time independent model around a steady state point for a model
global Tags modelParam WaitBar SimSettings IterCount TagInf TagFinal
Y0 =  modelParam.IC;
for n = 1:length(SetPoints)
    setPoint =SetPoints(n);
    %% find values at steady state point
%     SimSettings.PowerDemand = [SimSettings.PowerDemand(end),SimSettings.PowerDemand(end),setPoint,setPoint];
%     SimSettings.PowerTime = [0 .25 .5 1]*3600*24;
    SimSettings.PowerDemand = [SimSettings.PowerDemand(end),setPoint];
    SimSettings.PowerTime = [0 1]*3600*4;
    SimSettings.RunTime = SimSettings.PowerTime(end);

    IterCount = 1; TagInf =[]; TagFinal =[]; WaitBar.Show = 0;
%     WaitBar.Show = 1;  WaitBar.Text = 'Running non-linear model with transient';WaitBar.Handle =waitbar(0,WaitBar.Text);
    [~, Y] = ode15s(@RunBlocks, [0, SimSettings.RunTime],Y0);
%     close(WaitBar.Handle);
    Y0 = Y(end,:)';%states at steady state
    ControlStates = [];%find controller states
    controls = fieldnames(modelParam.Controls);
    for i = 1:length(controls)
        ControlStates((end+1):(end+length(modelParam.(controls{i}).States))) = modelParam.(controls{i}).States;
    end
    %% convert steady state inputs & outputs to vector and keep names of inputs/components
    U0 = InOutStates(Tags.ModelInput);
    Out0 = InOutStates(Tags.ModelOutput);

    %% find effect of state perturbations in the model portion
    ModelStates = (1:1:length(Y0))';
    ModelStates(ControlStates) = [];
    Y1 = Y0(ModelStates);
    Tags.Options.AssignedInputs =1;%start using assigned inputs in blocks instead of calculating new ones
    WaitBar.Show = 0;
    A = zeros(length(Y1));
    C = zeros(length(Out0));
    Perturbation = eye(length(Y0))*1e-9;
    for i = 1:length(Y1)%for all non-ignored states
        epsY = 2e-9;
        if ~(modelParam.Scale(ModelStates(i))==1 && Y1(i) ==1) % avoid valve that is full open. Perturbing in wrong direction will cause negatives
            dYplus = RunBlocks(SimSettings.RunTime,Y0+Perturbation(:,ModelStates(i)));%record resulting dY's
            dYplus(ControlStates)=[];%ignore dY's from ignored states
            OutPlus = InOutStates(Tags.ModelOutput);
        else
            dYplus = 0*Y1;
            OutPlus = Out0;
            epsY = 1e-9;
        end
        if ~(modelParam.Scale(ModelStates(i))==1 && Y1(i)<1e-9) % avoid valve that is full closed. Perturbing in wrong direction will cause negatives
            dYminus = RunBlocks(SimSettings.RunTime,Y0-Perturbation(:,ModelStates(i)));
            dYminus(ControlStates) = [];
            OutMinus = InOutStates(Tags.ModelOutput);
        else
            dYminus = 0*Y1;
            OutMinus = Out0;
            epsY = 1e-9;
        end
        A(:,i) = (dYplus-dYminus)/epsY;
        C(:,i) = (OutPlus-OutMinus)/epsY;
    end

    %% find effect of input perturbations on dY and Outputs
    B = zeros(length(Y1),length(U0));
    D = zeros(length(Out0),length(U0));
    Perturbation = eye(length(U0))*1e-9;
    for i = 1:length(U0)%for all inputs
        epsU = 2e-9;
        if U0(i)~=1 % avoid valve that is full open. Perturbing in wrong direction will cause negatives
            Tags.ModelInput = OverideInput(Tags.ModelInput,U0 + Perturbation(:,i));
            dYplus = RunBlocks(SimSettings.RunTime,Y0);%record resulting dY's
            dYplus(ControlStates)=[];%ignore dY's from ignored states
            OutPlus = InOutStates(Tags.ModelOutput);
        else
            dYplus = 0*Y1;
            OutPlus = Out0;
            epsU = 1e-9;
        end
        if U0(i)==0 % avoid valve that is full closed. Perturbing in wrong direction will cause negatives
            Tags.ModelInput = OverideInput(Tags.ModelInput,U0 - Perturbation(:,i));
            dYminus = RunBlocks(SimSettings.RunTime,Y0);
            dYminus(ControlStates) = [];
            OutMinus = InOutStates(Tags.ModelOutput);
        else
            dYminus = 0*Y1;
            OutMinus = Out0;
            epsU = 1e-9;
        end

        B(:,i) = (dYplus - dYminus)/epsU;
        D(:,i) = (OutPlus-OutMinus)/epsU;
    end

    m = round(length(SetPoints)/2);
    if n==m %create balanced realization and transform matrix from middle setpoint
        Sys = ss(A,B,C,D);% convert to model
        [Sys,g,T,Ti] = balreal(Sys);
        Model.Sys = Sys;
        Model.HSV = g;
        Model.T = T;
        Model.Ti = Ti;
    end
    
    Model.A = A;
    Model.B = B;
    Model.C = C;
    Model.D = D;

    Model.U0 = U0;
    Model.Out0 = Out0;
    Model.UX0 = Y0(ControlStates);
    Model.X = Y1;
    Model.GainScale = modelParam.Scale(ControlStates);
    Tags.Options = [];
    LinMod.Model{n} = Model;
    disp(strcat('Created linear model_',num2str(n),'_of_',num2str(length(SetPoints))))
end
LinMod = CatenateSys(LinMod,HSVtol);
controls = fieldnames(modelParam.Controls);
for i = 1:length(controls)
    LinMod.Controls.(controls{i}) = modelParam.(controls{i});
end
if isfield(modelParam,'Scope')
    LinMod.Scope = modelParam.Scope;
end
if isfield(modelParam,'NominalPower')
    LinMod.NominalPower = modelParam.NominalPower;
end
LinMod.IC = modelParam.IC;
    

function States = InOutStates(Input)
States = [];
list = fieldnames(Input);
for i = 1:1:length(list)
    block = list{i};
    ports = fieldnames(Input.(block));
    for j = 1:1:length(ports)
        port = ports{j};
        if isstruct(Input.(block).(port))
            f = fieldnames(Input.(block).(port));
            for k = 1:1:length(f)
                s = Input.(block).(port).(f{k});
                States(end+1:end+length(s),1) = s;
            end
        else
            s = Input.(block).(port);
            States(end+1:end+length(s),1) = s;
        end
    end
end

function Input = OverideInput(Input,New)
n=0;
list = fieldnames(Input);
for i = 1:1:length(list)
    block = list{i};
    ports = fieldnames(Input.(block));
    for j = 1:1:length(ports)
        port = ports{j};
        if isstruct(Input.(block).(port))
            f = fieldnames(Input.(block).(port));
            for k = 1:1:length(f)
                s = Input.(block).(port).(f{k});
                Input.(block).(port).(f{k}) = New(n+1:n+length(s));
                n = n+length(s);
            end
        else
            s = Input.(block).(port);
            Input.(block).(port) = New(n+1:n+length(s));
            n = n+length(s);
        end
    end
end

function Out = CatenateSys(System, Tol)
%removes states using balred and hankel singular values
l = length(System.Model);
m = round(l/2);
keep = (1:length(System.Model{1}.A));
keep(System.Model{m}.HSV < Tol) = [];
for i = 1:l
    A = System.Model{m}.T*System.Model{i}.A*System.Model{m}.Ti;
    B = System.Model{m}.T*System.Model{i}.B;
    C = System.Model{i}.C*System.Model{m}.Ti;
    X0 = System.Model{m}.T*System.Model{i}.X;

    Out.Model{i}.A = A(keep,keep);
    Out.Model{i}.B = B(keep,:);
    Out.Model{i}.C = C(:,keep);
    Out.Model{i}.D = System.Model{i}.D;
    Out.Model{i}.X0 = X0(keep);
    Out.Model{i}.Out0 = System.Model{i}.Out0;
    Out.Model{i}.U0 = System.Model{i}.U0;
    Out.Model{i}.UX0 = System.Model{i}.UX0;
    Out.Model{i}.GainScale = System.Model{m}.GainScale;
    Out.InterpVec(i) = X0(1);
end