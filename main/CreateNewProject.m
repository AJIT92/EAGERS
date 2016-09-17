function CreateNewProject(Name,hCreate)
global Plant SYSINDEX Model_dir figHandle Holidays
Plant.Name = Name;

%% Building Loads 
J = menu('Data Dialog','Generate Data From Building Models','Load Data From Files');
if J == 1
    h = msgbox('Automatically sets date for data to 2013');
    uiwait(h)
    Plant.Data.Timestamp = linspace(datenum(2013,1,1,.25,0,0),datenum(2014,1,1),96*(datenum(2014,1,1)-datenum(2013,1,1)));
    Plant.Data.Holidays = [735235;735255;735283;735318;735319;735320;735321;735322;735381;735402;735403;735404;735405;735406;735418;735479;735497;735498;735499;735500;735549;735566;735567;735591;735592;735593;735594;735599];%holidays of 2013, need to update
    Holidays = Plant.Data.Holidays;
    filename = fullfile(Model_dir, 'component library','Weather','1A.mat');
    load(filename)
    T=zeros(365*96,1);
    Temperature = (Temperature/10-32)*5/9; %convert temperature data to celcius
    T(1:4:end-3) = Temperature;
    T(2:4:end-2) = Temperature;
    T(3:4:end-1) = Temperature;
    T(4:4:end) = Temperature;
    Plant.Data.Temperature = T;
    Plant.Data.NameList = {};
    Plant.Data.Demand = {};
    Plant.Data.HistProf = {};
    Plant.Data.Weather = {};
    SetupBuilding
    waitfor(figHandle)
elseif J == 2
    LoadData()
    waitfor(figHandle)
end

Prompt = {'# of Electric Utilities','# of Gas Utilities','# of Electric Generators','#of CHP Generators','# of Chillers','# of Heaters','# of Boilers','# of Wind Generators', '# of Solar Generators', '# of Thermal Energy Storage Systems','# of Electric Energy Storage Systems'};
DefaultVal = {'1','1','1','0','0','0','0','0','0','0','0'};
A= inputdlg(Prompt,'Specify the micro-grid characteristics',1,DefaultVal);
SYSINDEX = 0;
Plant.Generator = [];
%% All Generators/heaters/chillers/boilers/wind/solar/
type = {'Electric Utility ';'Gas Utility';'Electric Generator';'CHP Generator';'Chiller';'Heater';'Boiler';'Wind';'Solar';'Thermal Storage';'Electric Storage';};%'Smart HVAC';};
for k = 1:1:length(type)
    for i = 1:1:str2double(A(k))
        SYSINDEX = SYSINDEX +1;
        if k ==1
            SetupUtility('Electric')
        elseif k==2
            SetupUtility('Gas')
        else
            J = menu(strcat('Specificify ',char(type(k)),' # ',num2str(i)),strcat('Load Existing ',char(type(k))),strcat('Modify Existing ',char(type(k))));
            MATfiles=dir(fullfile(Model_dir,'component library',char(type(k)),'*.mat'));
            list = {};
            for i = 1:1:length(MATfiles)
                list(i) = cellstr(MATfiles(i).name);
            end
            [s,v] = listdlg('PromptString','Select Component', 'SelectionMode','single','ListString',list);
            load(fullfile(Model_dir,'component library',char(type(k)),char(list(s))));
            Plant.Generator(SYSINDEX) =component;
            if J ==2
                if k<=7
                    SetupGenerator;
                elseif k==8
                    % SetupWind
                elseif k==9
                    SetupSolar   
                elseif k==10
                    SetupThermalStorage        
                elseif k==11
                    SetupBatteryStorage
                end
                waitfor(figHandle)
            end
        end
    end
end
%create options default
Plant.optimoptions.Interval = 1.0;% period of study in days
Plant.optimoptions.Horizon = 24;  %forecst window in hours
Plant.optimoptions.Resolution = 1.0; %interval frequency (may be a variable time step)
Plant.optimoptions.Topt = 600; % frequncy in s of the on-lne optimization
Plant.optimoptions.Tmpc = 60; % frequncy in s of the model predicitive controller loop
Plant.optimoptions.nsSmooth = 0; % temporal smoothing of heating demand
Plant.optimoptions.scaletime = 1.0; %time dilation so that 24 hour test can be caried out in shorter time (scales energy storage up)
Plant.optimoptions.fastsimulation = 1; % runs loops (dispatch, online, & MPC) sequentially rather than on timers
Plant.optimoptions.tspacing = 'constant';  %variable affecting time resolution, 1 means constant, 2 means linearly increasing, 3 means logarithmically increasing, 4 means manully input the values 
Plant.optimoptions.sequential = 1; %load all the chillers sequentially, so that they have their on fit curves for efficiency
Plant.optimoptions.excessHeat = 0; %excess heat from CCHP generators cannot be dumped to the environment. Heat demand is loaded into the equality equation.
Plant.optimoptions.thresholdSteps = 6;

%save plant
[f,p]=uiputfile(fullfile(Model_dir,'Plant',strcat(char(Plant.Name),'.mat')),'Save Project As...');
if f==0;return;end
save([p f],'Plant')
close(hCreate)