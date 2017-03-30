%% STRIDES %%
% this function can initialize and run a non-linear systems model
% It can also run a transient on that model
% It can also create a linear model, and compare the response for the same transient
global Model_dir SimSettings WaitBar modelParam LinMod IterCount TagInf TagFinal Tags
%% Load a Plant
%%%%-- either load a plant and initialize it, or load a saved Plant & model Parameters
%%%%-- Working options are: SOFCstack, SOECstack, SOFCsystem, GasTurbine
J = center_menu('Non-linear model options','Initialize from system description','Load previously initialized plant','Skip to pre-loaded linear model');
if J ==1
    %the goal is to have a GUI replace these m-files
    ModelFiles=dir(fullfile(Model_dir,'Model Library','*.m'));
    list = {};
    for i = 1:1:length(ModelFiles)
        list(i) = cellstr(strrep(ModelFiles(i).name,'.m',''));
    end
    [s,OK] = listdlg('PromptString','Select Model', 'SelectionMode','single','ListString',list);
    if OK~=1
        disp('Invalid selection. Exiting...')
    else
        modelName = list{s};
        Plant = feval(modelName);

        WaitBar.Show = 1;
        BuildModel(Plant); %%Build & Initialize model
        J2 = center_menu('Save Model?','Yes','No');
        if J2 ==1
            [f,p]=uiputfile(fullfile(Model_dir,'Model Library','Saved Models',strcat(modelName,'.mat')),'Save Model As...');
            save([p f],'modelParam')
        end
    end
elseif J ==2
    ModelFiles=dir(fullfile(Model_dir,'Model Library','Saved Models','*.mat'));
    list = {};
    for i = 1:1:length(ModelFiles)
        list(i) = cellstr(ModelFiles(i).name);
    end
    [s,OK] = listdlg('PromptString','Select Model', 'SelectionMode','single','ListString',list);
    if OK~=1
        disp('Invalid selection. Exiting...')
    else
        modelName = list{s};
        load(fullfile(Model_dir,'Model Library','Saved Models',modelName));
    end
end
%% Create or load a set of linear models
J2 = center_menu('Linear model option','Create new linearization','Load previously linearized plant','Skip linearization');
if J2==1 && J ==3
    disp('Error: You have not loaded a non-linear model to initialize')
    J2 =3;
end
if J2 ==1
    LinMod =[];
    HSVtol = 0; %what should this be set to?
    Prompt = {'Minimum Power (kW)','Maximum Power','# of Linear Models'};
    DefaultVal = {num2str(0.2*modelParam.NominalPower),num2str(modelParam.NominalPower),'5'};
    A= inputdlg(Prompt,'Specify range resolution of lineraization',1,DefaultVal);
    SetPoints = linspace(str2double(A(2)),str2double(A(1)),str2double(A(3)));
    LinMod = CreateLinModel(SetPoints,HSVtol);
    J3 = center_menu('Save Linearized Model?','Yes','No');
    if J3 ==1
        [f,p]=uiputfile(fullfile(Model_dir,'Model Library','Saved Linearizations',strcat(modelName,'.mat')),'Save Linearized Model As...');
        save([p f],'LinMod')
    end
elseif J2 ==2
    ModelFiles=dir(fullfile(Model_dir,'Model Library','Saved Linearizations','*.mat'));
    list = {};
    for i = 1:1:length(ModelFiles)
        list(i) = cellstr(ModelFiles(i).name);
    end
    [s,OK] = listdlg('PromptString','Select Model', 'SelectionMode','single','ListString',list);
    if OK~=1
        disp('Invalid selection. Exiting...')
    else
        load(fullfile(Model_dir,'Model Library','Saved Linearizations',list{s}));
        if isfield(LinMod,'NominalPower')
            modelParam.NominalPower = LinMod.NominalPower;
        else modelParam.NominalPower = 0;
        end
        modelParam.Controller = LinMod.Controls.Controller;
        modelParam.IC = LinMod.IC;
    end
elseif J2 ==3
    LinMod =[];
end

J3 = center_menu('Simulation Options','Simulate non-linear response','Simulate linear response','Simulate both linear and non-linear response','Neither');
if (J3 ==2 || J3 == 3) && J2 ==3
    disp('error you have not loaded a linearized model')
    J3 = 0;
end
if J3 ~=4 
    %set up the transient
    Prompt = {'Time','Demand (% of Nominal)'};
    DefaultVal = {'[0 4*3600 8*3600 24*3600]','[100 100 50 50]'};
    A= inputdlg(Prompt,'Specify the transient to test. Any string will be evaluated, but must create vertical vectors of equal length.',1,DefaultVal);
    SimSettings.PowerTime = eval(A{1});
    SimSettings.PowerDemand = eval(A{2})/100*modelParam.NominalPower;
    SimSettings.RunTime = SimSettings.PowerTime(end);
    

    %% modify control terms
    nC = length(modelParam.Controller.Gain);
    Prompt = {};
    DefaultVal = {};
    for i = 1:1:nC
        Prompt(end+1) = cellstr(strcat(modelParam.Controller.description{i},'---Proportional'));
        Prompt(end+1) = cellstr(strcat(modelParam.Controller.description{i},'---Integral'));
        DefaultVal(end+1) = cellstr(num2str(modelParam.Controller.PropGain(i)));
        DefaultVal(end+1) = cellstr(num2str(modelParam.Controller.Gain(i)));
    end
    A= inputdlg(Prompt,'Specify the controller parameters',1,DefaultVal);
    for i = 1:1:nC
        modelParam.Controller.PropGain(i) = str2double(A(2*i-1));
        modelParam.Controller.Gain(i) = str2double(A(2*i));
    end
    % %SOFCstack
    % modelParam.Controller.Gain = [3e-3;1e-3;1e-2];
    % modelParam.Controller.PropGain = [1;1;1];

    %SOFCsystem
    % modelParam.Controller.Gain = [1e-2;1e-4;1e-2];
    % modelParam.Controller.PropGain = [.5;.1;1];

    % %SOECstack
    % modelParam.Controller.Gain = [3e-3;1e-3;1e-2];
    % modelParam.Controller.PropGain = [1;1;1];

    % %GasTurbine
    % modelParam.Controller.IntGain = [4e-4; 1e-2; 4e-2;];
    % modelParam.Controller.PropGain = [8e-3; 5e-0; .75;];
end

%% Run a transient on non-linear model
if J3 ==1 || J3 == 3
    WaitBar.Show = 1;
    IterCount = 1; TagInf =[]; TagFinal =[];  WaitBar.Text = 'Running non-linear model with transient';WaitBar.Handle =waitbar(0,WaitBar.Text);
    tic; [T, Y] = ode15s(@RunBlocks, [0, SimSettings.RunTime], modelParam.IC); disp(strcat('Time to run model:',num2str(toc),' seconds'));close(WaitBar.Handle);
    PlotSimulation(T,Y,1,0,1)% Plot the tags and scopes in Plant.Plot. The first option (after Y) is to plot any fuel cell or electrolyzer temperature profiles, the second option is to mave a video of the transient, the third is to plot compressor turbine and blower maps
end

%% Run transient on set of linear models
if J3 ==2 || J3 ==3
    % need to find better initial condition (won't always start at nominal power)
    IC =  [LinMod.Model{1}.X0;LinMod.Model{1}.UX0];
%     Tags.O = LinMod.Model{1}.Out0;
    Tags.U = LinMod.Model{1}.U0;
%     Tags.dX = zeros(length(LinMod.Model{1}.X0),1);
%     Tags.dUX = zeros(length(LinMod.Model{1}.UX0),1);
    IterCount = 1; TagInf =[]; TagFinal =[]; WaitBar.Show = 1; WaitBar.Text = 'Running linear model with transient';WaitBar.Handle =waitbar(0,WaitBar.Text);
    tic; [T, Y] = ode15s(@RunLinSystem, [0, SimSettings.RunTime], IC); disp(strcat('Time to run model:',num2str(toc),' seconds'));close(WaitBar.Handle);
    if J3 == 3
        PlotSimulation(T,Y,0,0,0) % have already plotted the maps, just adding non-linear response to linear response
    else
        PlotSimulation(T,Y,1,0,1)
    end
end