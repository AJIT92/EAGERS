function SetupUtility(type)
global Plant SYSINDEX Model_dir figHandle
list = {};
if strcmp(type,'Electric')
    Plant.Generator(SYSINDEX).Type = 'Utility';
    K = menu('How would you like to specify the electric utility rate?','Peak/Off-Peak Rate','Load 15 Minute Data','Load Existing Utility','Modify Utility Structure');
    if K <=2
        Plant.Generator(SYSINDEX).Name = char(inputdlg('Electric Utility Name','Utility Name',1,{'Elec Utility1'}));
    end
    if K ==1
        if isfield (Plant.Generator(SYSINDEX), 'VariableStruct')
            if isfield(Plant.Generator(SYSINDEX).VariableStruct,'Erate')
            a = num2str(min(Plant.Generator(SYSINDEX).VariableStruct.Erate));
            b = num2str(max(Plant.Generator(SYSINDEX).VariableStruct.Erate));
            c = num2str(Plant.Generator(SYSINDEX).VariableStruct.MinThresh);
            d = num2str(max(Plant.Generator(SYSINDEX).VariableStruct.SellBack)/max(Plant.Generator(SYSINDEX).VariableStruct.Erate)*100);
            else a = '0.06';
            b = '0.08';
            c = '0';
            d = '0';
            end
            else a = '0.06';
            b = '0.08';
            c = '0';
            d = '0';
        end
        s = str2double(inputdlg({'Off-Peak electric rate in $ per kWh';'On-Peak (9AM-4PM) electric rate in $ per kWh';'Specify the minimum input threshold';'Specify the sellback rate (% of purchase rate)'},'Electric Rate', 1,{a;b;c;d;}));
        m = 96;% 15 min interval for utility pricing
        rate = [linspace(s(1),s(1),floor(8*m/24)) linspace(s(2),s(2),floor(12*m/24)) linspace(s(1),s(1),m-floor(20*m/24))]'; % $/kWh
        Erate = zeros(365*m+1,1); %full year
        for t = 1:1:length(Erate)
            Erate(t) = rate(1+mod(t-1,m));
        end
        MinThresh = s(3);
        SellBack = linspace(s(4),s(4),length(Erate))'.*Erate/100; % $/kWh
    elseif K == 2;
        %find xls files
        XLSfiles=dir(fullfile(Model_dir,'data','*.xls'));
        for i = 1:1:length(XLSfiles)
            list(i) = cellstr(XLSfiles(i).name);
        end
        %choose electric rate
        [s,v] = listdlg('PromptString','Select Electric Rate File (must be single column (15 min) on first sheet)', 'SelectionMode','single','ListString',list);
        Erate = xlsread(fullfile(Model_dir,'data',char(list(s)),'.xls'));
        s = str2double(inputdlg({'Specify the minimum input threshold';'Specify the sellback rate (% of purchase rate)'},'Electric Rate', 1,{'0';'0';}));
        MinThresh = s(1);
        SellBack = s(2)/100*Erate; % $/kWh
    end
    if K==1 || K==2
        Plant.Generator(SYSINDEX).Source = 'Electricity';
        Plant.Generator(SYSINDEX).Size = 1;
        Plant.Generator(SYSINDEX).VariableStruct.Erate = Erate;
        Plant.Generator(SYSINDEX).VariableStruct.MinImportThresh = MinThresh;
        Plant.Generator(SYSINDEX).VariableStruct.SellBack = SellBack;
        Plant.Generator(SYSINDEX).VariableStruct.WinStartMonth = 9;
        Plant.Generator(SYSINDEX).VariableStruct.WinStartDay = 22;
        Plant.Generator(SYSINDEX).VariableStruct.SumStartMonth = 3;
        Plant.Generator(SYSINDEX).VariableStruct.SumStartDay = 20;
        Plant.Generator(SYSINDEX).VariableStruct.SumRateTable = ones(7,24);
        Plant.Generator(SYSINDEX).VariableStruct.WinRateTable = ones(7,24);
        Plant.Generator(SYSINDEX).VariableStruct.SumRates = Erate;
        Plant.Generator(SYSINDEX).VariableStruct.WinRates = Erate;
    end
    if K==3 || K ==4
        %select a utility rate to modify
        MATfiles=dir(fullfile(Model_dir,'component library','Utility','*.mat'));
        for i = 1:1:length(MATfiles)
            list(i) = cellstr(MATfiles(i).name);
        end
        [s,v] = listdlg('PromptString','Select Electric Utility', 'SelectionMode','single','ListString',list);
        load(fullfile(Model_dir,'component library','Utility',char(list(s))));
        Plant.Generator(SYSINDEX).Name = component.Name;
        Plant.Generator(SYSINDEX).Source = 'Electricity';
        if K==3
             Plant.Generator(SYSINDEX).VariableStruct = component;
             figHandle = dialog('Visible','off');
             close(figHandle)
        elseif K ==4
            SetupElectricUtility(component,Plant.Data.Timestamp);
            waitfor(figHandle)
        end
    end
    Plant.Generator(SYSINDEX).Enabled = 1;
    Plant.Generator(SYSINDEX).Output.Capacity = linspace(0,1,11);
    Plant.Generator(SYSINDEX).Output.Electricity = linspace(0,1,11);
    Plant.Generator(SYSINDEX).Output.Heat = zeros(11,1);
    Plant.Generator(SYSINDEX).Output.Steam = zeros(11,1);
    Plant.Generator(SYSINDEX).Output.Cooling = zeros(11,1);
elseif strcmp(type,'Gas')
    Plant.Generator(SYSINDEX).Type = 'Utility';
    Plant.Generator(SYSINDEX).Name = char(inputdlg('Gas Utility Name','Utility Name',1,{'Gas Utility1'}));
    %% Specify gas rate
    K = menu('How would you like to specify the gas utility rate?','Constant Rate','Load 15 Minute Data');
    if K ==1
        if isfield(Plant.Generator(SYSINDEX).VariableStruct,'Grate')
            a= num2str(mean(Plant.Generator(SYSINDEX).VariableStruct.Grate));
        else a = num2str(4.15);
        end
        s = str2double(inputdlg('Specify the gas rate in $ per 1000 cu ft','Gas Rate', 1,{a}));
        a = 96*365+1;%1year @ 15 min resolution
        Grate = linspace(s,s,a)';
    elseif K == 2;
        %find xls files
        XLSfiles=dir(fullfile(Model_dir,'data','*.xls'));
        for i = 1:1:length(XLSfiles)
            list(i) = XLSfiles.name;
        end
        %choose electric rate
        [s,v] = listdlg('PromptString','Select Gas Rate File (must be single column (15 min) on first sheet)', 'SelectionMode','single','ListString',list);
        Grate = xlsread(fullfile(Model_dir,'data',list{s},'.xls'));
    end
    Plant.Generator(SYSINDEX).Source = 'NG';
    Plant.Generator(SYSINDEX).Output.Capacity = linspace(0,1,11);
    Plant.Generator(SYSINDEX).Output.Electricity = linspace(0,1,11);
    Plant.Generator(SYSINDEX).Output.Heat = linspace(0,1,11);
    Plant.Generator(SYSINDEX).Output.Steam = linspace(0,1,11);
    Plant.Generator(SYSINDEX).Output.Cooling = linspace(0,1,11);
    Plant.Generator(SYSINDEX).Size = 1;
    Plant.Generator(SYSINDEX).Enabled = 1;
    Plant.Generator(SYSINDEX).VariableStruct.Rate = Grate;
    Plant.Generator(SYSINDEX).VariableStruct.Timestamp = [0:length(Grate)-1].*(.25/24)+735235;
end