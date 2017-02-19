function varargout = NationalSurvey(varargin)
% NATIONALSURVEY M-file for NationalSurvey.fig
%      NATIONALSURVEY, by itself, creates a new NATIONALSURVEY or raises the existing
%      singleton*.
%
%      H = NATIONALSURVEY returns the handle to a new NATIONALSURVEY or the handle to
%      the existing singleton*.
%
%      NATIONALSURVEY('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in NATIONALSURVEY.M with the given input arguments.
%
%      NATIONALSURVEY('Property','Value',...) creates a new NATIONALSURVEY or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before NationalSurvey_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to NationalSurvey_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help NationalSurvey

% Last Modified by GUIDE v2.5 25-May-2014 23:30:45

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @NationalSurvey_OpeningFcn, ...
                   'gui_OutputFcn',  @NationalSurvey_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before NationalSurvey is made visible.
function NationalSurvey_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to NationalSurvey (see VARARGIN)

% Choose default command line output for NationalSurvey
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

global Project Model_dir
set(handles.uipanelSystem,'SelectedObject', handles.radiobuttonCurrentSys)
set(handles.uipanelUtility,'SelectedObject', handles.radiobuttonFixedRate)
set(handles.uipanelBuilding,'SelectedObject', handles.radiobuttonAllBuild)
set(handles.uipanelSysSize,'SelectedObject', handles.radiobuttonOptSize)
set(handles.editFileName,'string','Survey1')

compdir=fullfile(Model_dir, 'System Library', 'Solar');
files=dir(fullfile(compdir,'*.mat'));
list2=strrep({files.name},'.mat','');
list3(1) = {'Current Solar'};
for i = 1:1:length(list2)
    if (strcmp(list2(i),'DirectNormal')+strcmp(list2(i),'GlobalHorizontal')+strcmp(list2(i),'Azimuth')+strcmp(list2(i),'Zenith')) == 0
        list3(end+1) = list2(i);
    end
end
list2 = list3;
set(handles.popupmenuSolar,'string',list2)

compdir=fullfile(Model_dir, 'System Library', 'Wind');
files=dir(fullfile(compdir,'*.mat'));
list2=strrep({files.name},'.mat','');
list3(1) = {'Current Wind'};
for i = 1:1:length(list2)
    if (strcmp(list2(i),'WindPow')+strcmp(list2(i),'WindSpeed')) == 0
        list3(end+1) = list2(i);
    end
end
list2 = list3;
set(handles.popupmenuWind,'string',list2)

if isfield(Project,'Renewable')
    if isfield(Project.Renewable, 'Solar')
        AnnualPower = zeros(length(Project.Renewable.Solar),1);
        for i = 1:1:length(Project.Renewable.Solar)
            solar = Project.Renewable.Solar(i);
            if strcmp(solar.Tracking,'fixed')
                AnnualPower(i) = 8760*solar.Sizem2*sum((solar.Irrad/1000).*cosd(solar.SunZen-solar.Tilt).*(cosd(solar.SunAz-solar.Azimuth))*solar.Eff)/length(solar.SunAz);
            elseif strcmp(solar.Tracking,'1axis')
                AnnualPower(i) = 8760*solar.Sizem2*sum((solar.Irrad/1000).*cosd(solar.SunZen-solar.Tilt)*solar.Eff)/length(solar.SunAz);
            else AnnualPower(i) = 8760*solar.Sizem2*sum((solar.Irrad/1000)*solar.Eff)/length(solar.SunAz);
            end 
        end
        CHPSize = 0;
        for i = 1:1:length(Project.System.CHP)
            CHPSize = CHPSize+Project.System.CHP(i).SysSize(1);
        end
        Ts = 8760/length(Project.Building.DemandE);
        AnnualDem = sum(Project.Building.DemandE)*Ts;
        DemFrac = sum(AnnualPower)/AnnualDem;
        set(handles.editSizeAdem,'string',DemFrac*100)
    else 
        set(handles.editSizeAdem,'string',0)
    end
    if isfield(Project.Renewable, 'Wind')
        AnnualPower = zeros(length(Project.Renewable.Wind),1);
        for i = 1:1:length(Project.Renewable.Wind)
            wind = Project.Renewable.Wind(i);
            rho = 1.1798-1.3793e-4*wind.Elev+5.667e-9*wind.Elev^2;
            Wind = wind.Wind.*(wind.Wind>wind.CutIn).*(wind.Wind<wind.ShutDown);
            AnnualPower(i) = sum(wind.Eff*.5*rho*(3.1416*(wind.Diam)^2/4).*Wind.^3);
        end
        AnnualDem = sum(Project.Building.DemandE)*Ts;
        DemFrac = sum(AnnualPower)/AnnualDem;
        set(handles.editSizeAdem2,'string',DemFrac*100)
    else 
        set(handles.editSizeAdem2,'string',0)
    end
else
    set(handles.editSizeAdem,'string',0)
    set(handles.editSizeAdem2,'string',0)
end
set(handles.editSellBack,'string',num2str(Project.Utilities.Grid.SellBackRate))
if ~isfield(Project.Utilities.Grid,'SellBackPerc');
    Project.Utilities.Grid.SellBackPerc =100;
end
set(handles.editSellbackPerc,'string',num2str(Project.Utilities.Grid.SellBackPerc));
if Project.Utilities.Grid.SellBackRate == 0
    set(handles.uipanelSellBack,'SelectedObject',handles.radiobuttonNoSellBack)
elseif Project.Utilities.Grid.SellBackRate == -1
    set(handles.uipanelSellBack,'SelectedObject',handles.radiobuttonReverseMeter)
else
    set(handles.uipanelSellBack,'SelectedObject',handles.radiobuttonFixedSellBack)
end
set(handles.uipanelGridEmissions,'SelectedObject',handles.radiobuttonCombustionOnly)


% --- Outputs from this function are returned to the command line.
function varargout = NationalSurvey_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in popupmenuSolar.
function popupmenuSolar_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuSolar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Project
val=get(hObject,'value');
str=get(hObject,'string');
currentstr=str{val};
if val>1
    Project.Renewable = rmfield(Project.Renewable,'Solar');
    compdir=fullfile(Model_dir, 'System Library', 'Solar');
    load(fullfile(compdir,currentstr)) 
    Project.Renewable.Solar = component;
    AnnualPower = zeros(length(Project.Renewable.Solar),1);
    for i = 1:1:length(Project.Renewable.Solar)
        solar = Project.Renewable.Solar(i);
        if strcmp(solar.Tracking,'fixed')
            AnnualPower(i) = 8760*solar.Sizem2*sum((solar.Irrad/1000).*cosd(solar.SunZen-solar.Tilt).*(cosd(solar.SunAz-solar.Azimuth))*solar.Eff)/length(solar.SunAz);
        elseif strcmp(solar.Tracking,'1axis')
            AnnualPower(i) = 8760*solar.Sizem2*sum((solar.Irrad/1000).*cosd(solar.SunZen-solar.Tilt)*solar.Eff)/length(solar.SunAz);
        else AnnualPower(i) = 8760*solar.Sizem2*sum((solar.Irrad/1000)*solar.Eff)/length(solar.SunAz);
        end 
    end
    CHPSize = 0;
    for i = 1:1:length(Project.System.CHP)
        CHPSize = CHPSize+Project.System.CHP(i).SysSize(1);
    end

    AnnualDem = sum(Project.Building.DemandE)*Ts;
    DemFrac = sum(AnnualPower)/AnnualDem;
    set(handles.editSizeAdem,'string',DemFrac*100)
end

% --- Executes during object creation, after setting all properties.
function popupmenuSolar_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuSolar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenuWind.
function popupmenuWind_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuWind (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Project
val=get(hObject,'value');
str=get(hObject,'string');
currentstr=str{val};
if val>1
    Project.Renewable = rmfield(Project.Renewable,'Wind');
    compdir=fullfile(Model_dir, 'System Library', 'Wind');
    load(fullfile(compdir,currentstr)) 
    Project.Renewable.Wind = component;
    AnnualPower = zeros(length(Project.Renewable.Wind),1);
    for i = 1:1:length(Project.Renewable.Wind)
        wind = Project.Renewable.Wind(i);
        rho = 1.1798-1.3793e-4*wind.Elev+5.667e-9*wind.Elev^2;
        Wind = wind.Wind.*(wind.Wind>wind.CutIn).*(wind.Wind<wind.ShutDown);
        AnnualPower(i) = sum(wind.Eff*.5*rho*(3.1416*(wind.Diam)^2/4).*Wind.^3);
    end
    AnnualDem = sum(Project.Building.DemandE)*Ts;
    DemFrac = sum(AnnualPower)/AnnualDem;
    set(handles.editSizeAdem2,'string',DemFrac*100)
end

% --- Executes during object creation, after setting all properties.
function popupmenuWind_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuWind (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editSizeAdem_Callback(hObject, eventdata, handles)
% hObject    handle to editSizeAdem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Project
if isfield(Project.Renewable, 'Solar')
    SizeAdem = str2double(get(hObject,'String'));
    SolarSizing(1,'demand','all',SizeAdem)
    set(handles.editSizeAdem,'string',sum(Project.Renewable.Solar(:).DemFrac))
end

% --- Executes during object creation, after setting all properties.
function editSizeAdem_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editSizeAdem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function editSizeAdem2_Callback(hObject, eventdata, handles)
% hObject    handle to editSizeAdem2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Project
if isfield(Project.Renewable, 'Wind')
    SizeAdem2 = str2double(get(hObject,'String'));
    WindSizing(1,'demand','all',SizeAdem2)
    set(handles.editSizeAdem2,'string',sum(Project.Renewable.Wind(:).DemFrac))
end

% --- Executes during object creation, after setting all properties.
function editSizeAdem2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editSizeAdem2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pushbuttonCancel.
function pushbuttonCancel_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonCancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close(gcf)

function editFileName_Callback(hObject, eventdata, handles)
% hObject    handle to editFileName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object creation, after setting all properties.
function editFileName_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pushbuttonRun.
function pushbuttonRun_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonRun (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Project Model_dir 
h=waitbar(0,'Running Analysis');
Survey.StudyName = get(handles.editFileName,'String');
stateName = {'Alabama';'Alaska';'Arizona';'Arkansas';'California';'Colorado';'Connecticut';'Delaware';'Florida';'Georgia';
             'Hawaii';'Idaho';'Illinois';'Indiana';'Iowa';'Kansas';'Kentucky';'Louisiana';'Maine';'Maryland';
             'Massachusetts';'Michigan';'Minnesota';'Mississippi';'Missouri';'Montana';'Nebraska';'Nevada';'NewHampshire';'NewJersey';
             'NewMexico';'NewYork';'NorthCarolina';'NorthDakota';'Ohio';'Oklahoma';'Oregon';'Pennsylvania';'RhodeIsland';'SouthCarolina';
             'SouthDakota';'Tennessee';'Texas';'Utah';'Vermont';'Virginia';'Washington';'WestVirginia';'Wisconsin';'Wyoming';};
stateAbrev = {'AL';'AK';'AZ';'AR';'CA';'CO';'CT';'DE';'FL';'GA';'HI';'ID';'IL';'IN';'IA';'KS';'KY';'LA';'ME';'MD';'MA';'MI';'MN';'MS';'MO';
              'MT';'NE';'NV';'NH';'NJ';'NM';'NY';'NC';'ND';'OH';'OK';'OR';'PA';'RI';'SC';'SD';'TN';'TX';'UT';'VT';'VA';'WA';'WV';'WI';'WY';}; 
buildTypeName = {'SDRest'; 'FFRest'; 'Sch-pri'; 'Sch-sec'; 'LgOff'; 'MdOff'; 'SmOff'; 'MRapt'; 'LgHotel'; 'SmHotel'; 'Hospital'; 'OutP'; 'Retail'; 'StMall'; 'SMarket'; 'ware';};
 
if get(handles.uipanelSystem,'SelectedObject')~= handles.radiobuttonCurrentSys
    load(fullfile(fullfile(Model_dir,'project'),'DefaultSystem.mat'))
    %% load analysis type (set System, FixSize, Control & heat recovery)
    if get(handles.uipanelSystem,'SelectedObject') == handles.radiobuttonFC
        if isfield(Project.System,'Chiller')
            A = rmfield(Project.System,'Chiller');
            Project.System = A;
        end
        if isfield(Project.System,'TES')
            A = rmfield(Project.System,'TES');
            Project.System = A;
        end
    elseif get(handles.uipanelSystem,'SelectedObject') == handles.radiobuttonFCchill
        A = rmfield(Project.System,'TES');
        Project.System = A;
    end
end
if isfield(Project,'Renewable')
    if isfield(Project.Renewable,'Solar')
        SolarDir = fullfile(Model_dir,'System Library' ,'Solar'); 
        load(strcat(SolarDir,[filesep 'solarData' filesep 'GlobalHorizontal.mat']));
        load(strcat(SolarDir,[filesep 'solarData' filesep 'DirectNormal.mat']));
        load(strcat(SolarDir,[filesep 'solarData' filesep 'Azimuth.mat']));
        load(strcat(SolarDir,[filesep 'solarData' filesep 'Zenith.mat']));
    end
    if isfield(Project.Renewable,'Wind')
        WindDir = fullfile(Model_dir,'System Library', 'Wind');
        load(strcat(WindDir,[filesep 'windData' filesep 'WindSpeed.mat']));
        load(strcat(WindDir,[filesep 'windData' filesep 'WindPow.mat']));
    end
end

BuildDir = fullfile(Model_dir,'System Library','Buildings');
building = Project.Building;
if get(handles.uipanelBuilding,'SelectedObject')== handles.radiobuttonAllBuild
    Survey.BuildName = 'all';
else
    for i = 1:1:length(building.NameList)
        name = char(building.NameList{i});
        ind=strfind(name,'_');
        Survey.BuildName(i) = cellstr(name(1:ind(1)-1));
    end
end

A = get(handles.uipanelGridEmissions,'SelectedObject');
switch get(A,'Tag')
    case 'radiobuttonStateAvg'
        GridMix = 1;
    case 'radiobuttonCombustionOnly'
        GridMix = 2;
end
if get(handles.uipanelSysSize,'SelectedObject') == handles.radiobuttonFixSize
    SizeStrategy = 1; %dont resize
elseif get(handles.uipanelSysSize,'SelectedObject') == handles.radiobuttonOptSize
    SizeStrategy = 2;
elseif get(handles.uipanelSysSize,'SelectedObject') == handles.radiobuttonBestNPV
    SizeStrategy = 3;
elseif get(handles.uipanelSysSize,'SelectedObject') == handles.radiobuttonLowestCO2
    SizeStrategy = 4;
end
BaseSystem = Project.System;
if isfield(Project,'Renewable')
    BaseRenewables = Project.Renewable;
else BaseRenewables =[];
end
for state = 1:1:50
    %% Change State
    Project.State = stateName(state);
    if ~isempty(BaseRenewables)
        Project.Renewable = BaseRenewables;
    end
    %%load new building profiles (Change climate)
    if get(handles.uipanelBuilding,'SelectedObject')== handles.radiobuttonAllBuild
        buildTest = length(buildTypeName);
    else buildTest =1;
    end
    for j = 1:1:buildTest
        if get(handles.uipanelBuilding,'SelectedObject')== handles.radiobuttonAllBuild
            name = char(building.Name);
            if strcmp(name, 'multiple')
                name = strcat(buildTypeName(1),'_',char(ClimateZone(stateAbrev(state),0 )),'_','New2010');
            end
            ind=strfind(name,'_');
            oldBuildType = name(1:ind(1)-1);
            OldClimate = name(ind(1):ind(2));
            name = strrep(name,oldBuildType,char(buildTypeName(j)));
            name = strrep(name,OldClimate,strcat('_',char(ClimateZone(stateAbrev(state),0 )),'_'));
            load(strcat(BuildDir,name))
            %%% edit exterior lighting profile
            AbsMin =min(component.DemandE);
            LightLoad = max(component.ExteriorLight);
            AddLight = (component.ExteriorLight==0).*(component.DemandE<=(AbsMin+LightLoad));
            component.DemandE = component.DemandE+AddLight*LightLoad;
            %%% 
            building.Name = name;
            building.NameList(1) = cellstr(name);
            building.DemandE = component.DemandE;
            building.DemandC = building.DistCool*component.DemandC;
            building.DemandH = building.DistHeat*component.DemandH;
            building.NonDistH =(1-building.DistHeat)*component.DemandH;
            building.CoolingElectricalLoad =  building.DistCool*component.CoolingElectricalLoad;
        elseif get(handles.uipanelBuilding,'SelectedObject')== handles.radiobuttonCurrent
            building.DemandE = building.DemandE*0;
            building.DemandC = building.DemandC*0;
            building.DemandH = building.DemandH*0;
            building.CoolingElectricalLoad = building.DemandE*0;
            for i = 1:1:length(building.NameList)
                name = char(building.NameList{i});
                ind=strfind(name,'_');
                OldClimate = name(ind(1):ind(2));
                name = strrep(name,OldClimate,strcat('_',char(ClimateZone(stateAbrev(state),0)),'_'));
                load(strcat(BuildDir,name))
                %%% edit exterior lighting profile
                AbsMin =min(component.DemandE);
                LightLoad = max(component.ExteriorLight);
                AddLight = (component.ExteriorLight==0).*(component.DemandE<=(AbsMin+LightLoad));
                component.DemandE = component.DemandE+AddLight*LightLoad;
                %%% 
                building.DemandE = building.DemandE + component.DemandE;
                building.DemandC = building.DemandC + building.DistCool(i)*component.DemandC;
                building.DemandH = building.DemandH + building.DistHeat(i)*component.DemandH;
                building.NonDistH = building.NonDistH + (1-building.DistHeat(i))*component.DemandH;
                building.CoolingElectricalLoad = building.CoolingElectricalLoad + building.DistCool(i)*component.CoolingElectricalLoad;
                building.NameList(i) = cellstr(name);
            end
            if length(building.NameList) ==1
                building.Name = name;
            else building.Name = 'multiple';
            end
        end
    Project.Building = building;
    %End Loading Building
    %% load new renewables
    if j==1 && isfield(Project,'Renewable')
        
        if isfield(Project.Renewable,'Solar')
            for i=1:1:length(Project.Renewable.Solar)
                Project.Renewable.Solar(i).State=stateName(state);
                if strcmp(Project.Renewable.Solar(i).PVtype,'flat')
                    Project.Renewable.Solar(i).Irrad = GlobalHorizontal(:,state);
                else Project.Renewable.Solar(i).Irrad = DirectNormal(:,state);
                end
                Project.Renewable.Solar(i).SunAz = Azimuth(:,state);
                Project.Renewable.Solar(i).SunZen = Zenith(:,state);
                if strcmp(Project.Renewable.Solar(i).Tracking,'fixed')
                    track = 1;
                elseif strcmp(Project.Renewable.Solar(i).Tracking,'1axis')
                    track = 2;
                else track =3;
                end
                Project.Renewable.Solar(i).Tilt = OptimalZenith(Project.Renewable.Solar(i).Irrad,Project.Renewable.Solar(i).SunAz,Project.Renewable.Solar(i).SunZen,track,Project.Renewable.Solar(i).Azimuth,Project.Renewable.Solar(i).Eff);
            end
        end
        if isfield(Project.Renewable,'Wind')
            for i=1:1:length(Project.Renewable.Wind)
                Project.Renewable.Wind(i).State=stateName(state);
                Project.Renewable.Wind(i).Wind = WindSpeed(:,state);
            end
        end
    end
    Project.System = BaseSystem;
    %% Re-size applicable parts of system
    if  SizeStrategy==1 % Resize system components if sizing strategy selected
        %do nothing, fixed size
    elseif SizeStrategy==2
        ScaleSystemComponents(100,100,100,120,1,2,1,2,2)
    elseif SizeStrategy==3
        if isfield(Project.System,'TES')
            autoSizing('cost','TES')
        end
        autoSizing('cost','CHP')
        if isfield(Project.System,'Battery')
            autoSizing('cost','Battery',3)
        end
    elseif SizeStrategy==4
        autoSizing('CO2','CHP')
    end
    %% Re-size applicable parts of system
    if get(handles.uipanelUtility,'SelectedObject') == handles.radiobuttonFixedRate
        scale=0;
    else scale =1;
    end
    %% Run Project
    Project.Result = runAnalyses(Project,'Commercial',scale);

    %for use in results fig
    Survey.Econ = Project.Economic;
    Survey.Sizing(state,j).Sizes = [Project.Result.Dispatch.SysSize Project.Result.Dispatch.ChillerSize(1) Project.Result.Dispatch.ChillerSize(2) Project.Result.Dispatch.TESsize Project.Result.Dispatch.BatterySize Project.Result.Dispatch.RenewableSize(1) Project.Result.Dispatch.RenewableSize(2)];
    Cost = Project.Result.costOut;
    %Determine break-even cost of FC system
    Survey.FCevenCost(state,j) = BreakEvenFCcost(Project.Economic,Cost,Project.Result.eOut.SysSize);
    Survey.FCsize(state,j) = Project.Result.eOut.SysSize;
    %% record other variables 
    Survey.NPCdispatch(1:5,state,j) = [Cost.NPVnewDemCharges Cost.NPVnewUseCharges Cost.NPVnewFuelCost Cost.NPVnewOandM Cost.NPVnewFinance];
    Survey.NPCbaseline(1:5,state,j) = [Cost.NPVbaselineDemCharges Cost.NPVbaselineUseCharges Cost.NPVbaselineFuelCost Cost.NPVbaselineOandM Cost.NPVbaselineFinance];
    Survey.Payback(state,j) = Cost.Payback;
    Survey.Year1baseline(1:3,state,j) = Cost.Year1baseCharges;
    Survey.Year1dispatch(1:3,state,j) = Cost.Year1dispatchCharges;
    Survey.SelfGen(state,j) = Project.Result.eOut.SelfGen;
    Survey.CapacitykW(state,j) = Project.Result.Dispatch.SysSize;
    Survey.GenkWh(state,j) = Project.Result.Dispatch.ElecTotProd;
    Survey.CapacityFactor(state,j) = Survey.GenkWh(state,j)/(8760*Survey.CapacitykW(state,j));

    % Emissions
    [BaselineEmission DispatchEmission] = EmissionsCalculated(Project.Result,stateAbrev(state),GridMix,1);

    Survey.Dispatch(state,j).CO2 = sum(DispatchEmission(1,:));
    Survey.Dispatch(state,j).NOx = sum(DispatchEmission(2,:));
    Survey.Dispatch(state,j).SO2 = sum(DispatchEmission(3,:));
    Survey.Baseline(state,j).CO2= sum(BaselineEmission(1,:));
    Survey.Baseline(state,j).NOx = sum(BaselineEmission(2,:));
    Survey.Baseline(state,j).SO2 = sum(BaselineEmission(3,:));
    
    % Demand Characteristics (describe shape of demand profile
    Survey.LoadFactor(state,j) = mean(Project.Building.DemandE)/max(Project.Building.DemandE);
    Survey.LoadScatter(state,j) = std(Project.Building.DemandE)/mean(Project.Building.DemandE);
    Survey.LoadEtoC(state,j) = mean(Project.Building.DemandE)/mean(Project.Building.DemandC);
    Survey.LoadEtoH(state,j) = mean(Project.Building.DemandE)/mean(Project.Building.DemandH);
    Survey.LoadCtoH(state,j) = mean(Project.Building.DemandC)/mean(Project.Building.DemandH);

    waitbar(state/50,h,strcat('Running Analysis',num2str(j+buildTest*(state-1)),' of ',num2str(50*buildTest)));
    end
end
close(h)

ResultDir=fullfile(Model_dir, 'results'); 
name = strcat('NatSurv',get(handles.editFileName,'String'));
filename=fullfile(ResultDir,name);
save(filename,'Survey')
close(gcf)
%% Visualization                
NationalSurveyResults()

function Tilt = OptimalZenith(Irrad,Az,Zen,track,Ao,Eff)
minZen = min(Zen+99*(Zen==0));
meanZen = sum(Zen)/sum(Zen>0);

% Test multiple tilt angles
To = linspace(minZen,meanZen,20);
Col = linspace(1,1,length(To));
if track ==3
    Tilt = 0;
elseif track ==1
    Po = (Irrad*Col).*cosd(Zen*Col-linspace(1,1,length(Zen))'*To).*(cosd(Az-Ao)*Col)*Eff;
    AnPow = sum(Po,1);
    [MaxPow, I] = max(AnPow);
    if I ==1
        I =2;
    end
    To = linspace(To(I-1),To(I+1),20);
    Po = (Irrad*Col).*cosd(Zen*Col-linspace(1,1,length(Zen))'*To).*(cosd(Az-Ao)*Col)*Eff;
    AnPow = sum(Po,1);
    [MaxPow, I] = max(AnPow);
    Tilt = To(I);
elseif track ==2
    Po = (Irrad*Col).*cosd(Zen*Col-linspace(1,1,length(Zen))'*To)*Eff;
    AnPow = sum(Po,1);
    [MaxPow, I] = max(AnPow);
    if I ==1
        I =2;
    end
    To = linspace(To(I-1),To(I+1),20);
    Po = (Irrad*Col).*cosd(Zen*Col-linspace(1,1,length(Zen))'*To).*(cosd(Az-Ao)*Col)*Eff;
    AnPow = sum(Po,1);
    [MaxPow, I] = max(AnPow);
    Tilt = To(I);
end

function basePerc = baseLoad(demand,grid)
steps = length(demand.DemandE);
Ts = 8760/steps; 
%Find holidays/weekends
ElecOnly = demand.DemandE-demand.CoolingElectricalLoad;
dayMax = zeros(365,1);
for day = 1:1:365
    dayMax(day) = max(ElecOnly(1+24/Ts*(day-1):24/Ts*day));
end
holidays = [];
for i = 1:52
    lowDays =7*(i-1) + find(dayMax(1+7*(i-1):7*i)<(median(dayMax(1+7*(i-1):7*i))-.5*std(dayMax(1+7*(i-1):7*i))));
    holidays(end+1:end+length(lowDays)) = lowDays;
end

%Find Minimum demand without weekends/Holidays
SumStart = datenum(2013,grid.summerStartDate(1), grid.summerStartDate(2))-datenum(2013,1,1);
SumEnd = datenum(2013,grid.summerEndDate(1), grid.summerEndDate(2))-datenum(2013,1,1);
SummerHolidays = nnz((holidays>=SumStart).*(holidays<=SumEnd));
count = 0;
dayProf= zeros((SumEnd-SumStart+1)-SummerHolidays,24/Ts);
for day = SumStart:1:SumEnd  %%Size system for summer loads
    if nnz(day == holidays)==0;
        count = count+1;
        dayProf(count,1:24/Ts) = demand.DemandE(1+24/Ts*(day-1):24/Ts*day);
    end
end
MinDemandE = min(min(dayProf));
basePerc = 8760*MinDemandE/sum(demand.DemandE*Ts)*100;




function editSellBack_Callback(hObject, eventdata, handles)
% hObject    handle to editSellBack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Project.Utilities.Grid.SellBackRate = str2double(get(hObject,'String'));
set(handles.uipanelSellBack,'SelectedObject',handles.radiobuttonFixedSellBack)

% --- Executes during object creation, after setting all properties.
function editSellBack_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editSellBack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editSellbackPerc_Callback(hObject, eventdata, handles)
% hObject    handle to editSellbackPerc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Project
Project.Utilities.Grid.SellBackPerc =str2double(get(hObject,'string'));


% --- Executes during object creation, after setting all properties.
function editSellbackPerc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editSellbackPerc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes when selected object is changed in uipanelChillType.
function uipanelSellBack_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanelChillType 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
global Project 
switch get(eventdata.NewValue,'Tag')
    case 'radiobuttonNoSellBack'
        Project.Utilities.Grid.SellBackRate = 0;
        set(handles.editSellBack,'string','0.00')
    case 'radiobuttonFixedSellBack'
        if Project.Utilities.Grid.SellBackRate < 0;
            Project.Utilities.Grid.SellBackRate = 0;
        else
            Project.Utilities.Grid.SellBackRate = str2double(get(handles.editSellBack,'String'));
        end
        set(handles.editSellBack,'string',num2str(Project.Utilities.Grid.SellBackRate))
    case 'radiobuttonReverseMeter'
        Project.Utilities.Grid.SellBackRate = -1;
        set(handles.editSellBack,'string','-1')
end


function FCcost = BreakEvenFCcost(Econ,Result,SysSize)
InterestMonth = Econ.Interest/12/100;
months = Econ.FinanceYrs*12;

BaseNPC = Result.NPVbaselineUseCharges + Result.NPVbaselineDemCharges + Result.NPVbaselineFuelCost + Result.NPVbaselineOandM + Result.NPVbaselineFinance;
Non_FC_NPC = Result.NPVnewUseCharges + Result.NPVnewDemCharges + Result.NPVnewFuelCost + Result.NPVnewOandM + Result.NonFCfinanceNPV;
NPCforFC =  BaseNPC-Non_FC_NPC;

Yearl1Cost = reverseNPC(NPCforFC,1+Econ.Inflation/100,Econ.FinanceYrs);
SystemCost = Yearl1Cost*(1-(1+InterestMonth)^(-months))/(12*InterestMonth);
FCcost = (SystemCost)/SysSize + Econ.Incentive;


function Yearl1Cost = reverseNPC(NPC,irr,years)
NPCtrend = zeros(years,1);
NPCtrend(1)=1/(irr);
for i=2:years
    NPCtrend(i)=NPCtrend(i-1)+1/(irr)^i;
end
Yearl1Cost = NPC/NPCtrend(end);


function value=NPC(irr,cashflows)

value=0;
for i=1:length(cashflows)
    value=value+cashflows(i)/(irr)^i;
end