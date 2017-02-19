function varargout = FleetAnalysis(varargin)
% FLEETANALYSIS M-file for FleetAnalysis.fig
%      FLEETANALYSIS, by itself, creates a new FLEETANALYSIS or raises the existing
%      singleton*.
%
%      H = FLEETANALYSIS returns the handle to a new FLEETANALYSIS or the handle to
%      the existing singleton*.
%
%      FLEETANALYSIS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FLEETANALYSIS.M with the given input arguments.
%
%      FLEETANALYSIS('Property','Value',...) creates a new FLEETANALYSIS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before FleetAnalysis_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to FleetAnalysis_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help FleetAnalysis

% Last Modified by GUIDE v2.5 21-Jan-2014 12:49:22

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @FleetAnalysis_OpeningFcn, ...
                   'gui_OutputFcn',  @FleetAnalysis_OutputFcn, ...
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


% --- Executes just before FleetAnalysis is made visible.
function FleetAnalysis_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to FleetAnalysis (see VARARGIN)

% Choose default command line output for FleetAnalysis
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

global Model_dir Fleet
%change to find project directory more reliably
fleetdir=fullfile(Model_dir, 'System Library', 'Fleet');
files=dir(fullfile(fleetdir,'*.mat'));
list=strrep({files.name},'.mat','');
load(fullfile(fleetdir,list{1}))
% load(fullfile(fleetdir,'DefaultFleet'))
set(handles.popupmenuFleet,'string',list,'value',1)

set(handles.uipanelSystem,'SelectedObject', handles.radiobuttonFC)
set(handles.uipanelControl,'SelectedObject', handles.radiobuttonBaseload)
set(handles.uipanelCHP,'SelectedObject', handles.radiobuttonHeatRecov)


set(handles.editSaveName,'string','WebStudy1')
set(handles.editBuildName,'string',Fleet.BuildList(1))
set(handles.editAvgLoad,'string',num2str(Fleet.Size(1)))

stateName = {'Alabama';'Alaska';'Arizona';'Arkansas';'California';'Colorado';'Connecticut';'Delaware';'Florida';'Georgia';
             'Hawaii';'Idaho';'Illinois';'Indiana';'Iowa';'Kansas';'Kentucky';'Louisiana';'Maine';'Maryland';
             'Massachusetts';'Michigan';'Minnesota';'Mississippi';'Missouri';'Montana';'Nebraska';'Nevada';'NewHampshire';'NewJersey';
             'NewMexico';'NewYork';'NorthCarolina';'NorthDakota';'Ohio';'Oklahoma';'Oregon';'Pennsylvania';'RhodeIsland';'SouthCarolina';
             'SouthDakota';'Tennessee';'Texas';'Utah';'Vermont';'Virginia';'Washington';'WestVirginia';'Wisconsin';'Wyoming';};
stateAbrev = {'AL';'AK';'AZ';'AR';'CA';'CO';'CT';'DE';'FL';'GA';'HI';'ID';'IL';'IN';'IA';'KS';'KY';'LA';'ME';'MD';'MA';'MI';'MN';'MS';'MO';
              'MT';'NE';'NV';'NH';'NJ';'NM';'NY';'NC';'ND';'OH';'OK';'OR';'PA';'RI';'SC';'SD';'TN';'TX';'UT';'VT';'VA';'WA';'WV';'WI';'WY';};   
BuildType = {'Restaurant: full-service (sit down)';'Restaurant: quick-service (fast food)';'School: primary school';'School: secondary school';'Office: large office';'Office: medium office';'Office: small office';'Mid-rise apartment building';'Hospitality: large hotel';'Hospitality: small hotel/motel';'Health care: large hospital';'Health care: outpatient facility';'Retail: big-box, standalone retail store';'Retail: retail store located in a strip mall';'Retail: supermarket';'Unrefrigerated warehouse';};
buildTypeName = {'SDRest'; 'FFRest'; 'Sch-pri'; 'Sch-sec'; 'LgOff'; 'MdOff'; 'SmOff'; 'MRapt'; 'LgHotel'; 'SmHotel'; 'Hospital'; 'OutP'; 'Retail'; 'StMall'; 'SMarket'; 'ware';};
bT = find(strcmp(Fleet.BuildType(1),buildTypeName));
sN = find(strcmp(Fleet.State(1),stateAbrev));
set(handles.lbFleet, 'string', Fleet.BuildList(:,1), 'Max', length(Fleet.BuildList(:,1)), 'Min', 1)
set(handles.popupmenuBuildType, 'string', BuildType(:,1), 'value',bT)
set(handles.popupmenuState, 'string', stateName(:,1), 'value',sN)
set(handles.editZipCode, 'string', '00000')
set(handles.editSolar, 'string', num2str(Fleet.SolarSize(1)))
lbFleet_Callback(handles.lbFleet, eventdata, handles)
      

% --- Outputs from this function are returned to the command line.
function varargout = FleetAnalysis_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
varargout{1} = handles.output;

% --- Executes on button press in pushbuttonCancel.
function pushbuttonCancel_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonCancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close(gcf)

% --- Executes on button press in pushbuttonFleet.
function pushbuttonFleet_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonFleet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Project Model_dir
projdir=fullfile(Model_dir,'System Library','Fleet');
[f,p]=uiputfile('*.mat','Save Fleet As...',fullfile(projdir, 'Fleet.mat'));
if f==0; return; end
Fleet=Project.Fleet;
save([p,f],'Fleet')

% --- Executes on selection change in popupmenuFleet.
function popupmenuFleet_Callback(hObject, eventdata, handles, skipreload)
% hObject    handle to popupmenuFleet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Project Model_dir
FleetList=get(hObject,'string');
FleetVal=get(hObject,'value');
if nargin ==3 || skipreload == 0
    fleetdir=fullfile(Model_dir, 'System Library','Fleet');
    if strcmp(FleetList{FleetVal},'DefaultBuilding')
        load(fullfile(fleetdir,FleetList{FleetVal}));
    else
        fleetName=fullfile(fleetdir,FleetList{FleetVal});
        load(fleetName)
        Project.Fleet = Fleet;
    end
end
set(handles.lbFleet,'string',Fleet.BuildList,'value',1)

% --- Executes on selection change in lbFleet.
function lbFleet_Callback(hObject, eventdata, handles)
% hObject    handle to lbFleet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Fleet
val=get(hObject,'value');
stateName = {'Alabama';'Alaska';'Arizona';'Arkansas';'California';'Colorado';'Connecticut';'Delaware';'Florida';'Georgia';
             'Hawaii';'Idaho';'Illinois';'Indiana';'Iowa';'Kansas';'Kentucky';'Louisiana';'Maine';'Maryland';
             'Massachusetts';'Michigan';'Minnesota';'Mississippi';'Missouri';'Montana';'Nebraska';'Nevada';'NewHampshire';'NewJersey';
             'NewMexico';'NewYork';'NorthCarolina';'NorthDakota';'Ohio';'Oklahoma';'Oregon';'Pennsylvania';'RhodeIsland';'SouthCarolina';
             'SouthDakota';'Tennessee';'Texas';'Utah';'Vermont';'Virginia';'Washington';'WestVirginia';'Wisconsin';'Wyoming';};
stateAbrev = {'AL';'AK';'AZ';'AR';'CA';'CO';'CT';'DE';'FL';'GA';'HI';'ID';'IL';'IN';'IA';'KS';'KY';'LA';'ME';'MD';'MA';'MI';'MN';'MS';'MO';
              'MT';'NE';'NV';'NH';'NJ';'NM';'NY';'NC';'ND';'OH';'OK';'OR';'PA';'RI';'SC';'SD';'TN';'TX';'UT';'VT';'VA';'WA';'WV';'WI';'WY';};
buildTypeName = {'SDRest'; 'FFRest'; 'Sch-pri'; 'Sch-sec'; 'LgOff'; 'MdOff'; 'SmOff'; 'MRapt'; 'LgHotel'; 'SmHotel'; 'Hospital'; 'OutP'; 'Retail'; 'StMall'; 'SMarket'; 'ware';};
BuildType = {'Restaurant: full-service (sit down)';'Restaurant: quick-service (fast food)';'School: primary school';'School: secondary school';'Office: large office';'Office: medium office';'Office: small office';'Mid-rise apartment building';'Hospitality: large hotel';'Hospitality: small hotel/motel';'Health care: large hospital';'Health care: outpatient facility';'Retail: big-box, standalone retail store';'Retail: retail store located in a strip mall';'Retail: supermarket';'Unrefrigerated warehouse';};
climateName = {'_1A_'; '_2A_'; '_2B_'; '_3A_'; '_3B_'; '_3B-Coast_'; '_3C_'; '_4A_'; '_4B_'; '_4C_'; '_5A_'; '_5B_'; '_6A_'; '_6B_'; '_7_'; '_8_';};
Climate = {'Miami (ASHRAE 1A)';'Houston (ASHRAE 2A)';'Phoenix (ASHRAE 2B)';'Atlanta (ASHRAE 3A)';'Las Vegas (ASHRAE 3B-Inland)';'Los Angeles (ASHRAE 3B-Coast)';'San Francisco (ASHRAE 3C)';'Baltimore (ASHRAE 4A)';'Albuquerque (ASHRAE 4B)';'Seattle (ASHRAE 4C)';'Chicago (ASHRAE 5A)';'Boulder (ASHRAE 5B)';'Minneapolis (ASHRAE 6A)';'Helena, MT (ASHRAE 6B)';'Duluth, MN (ASHRAE 7)';'Fairbanks, AK (ASHRAE 8)';};


I = find(strcmp(Fleet.BuildType(val),buildTypeName));
ZipCode = Fleet.ZipCode(val);
[ASHRAEZone] = ClimateZone( char(Fleet.State(val)),ZipCode );
K = find(strcmp(strcat('_',ASHRAEZone,'_'),climateName));
State = find(strcmp(Fleet.State(val),stateAbrev));
description = {strcat('Building Type:  ',char(BuildType(I,1)));
                strcat('Climate Type:  ',char(Climate(K,1)));
                strcat('State:  ',char(stateName(State)))
                strcat('Average Electric Demand (kW):  ',num2str(Fleet.Size(val,1)));
                strcat('Rooftop Photovoltaics (m2):  ',num2str(Fleet.SolarSize(val)));};
set(handles.textBuildDescrip,'string',description)

% --- Executes during object creation, after setting all properties.
function lbFleet_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lbFleet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pushbuttonAdd.
function pushbuttonAdd_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonAdd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Fleet
nFleetBuild = length(Fleet.BuildList);
type=get(handles.popupmenuBuildType, 'value');
State=get(handles.popupmenuState, 'value');
buildTypeName = {'SDRest'; 'FFRest'; 'Sch-pri'; 'Sch-sec'; 'LgOff'; 'MdOff'; 'SmOff'; 'MRapt'; 'LgHotel'; 'SmHotel'; 'Hospital'; 'OutP'; 'Retail'; 'StMall'; 'SMarket'; 'ware';};
stateAbrev = {'AL';'AK';'AZ';'AR';'CA';'CO';'CT';'DE';'FL';'GA';'HI';'ID';'IL';'IN';'IA';'KS';'KY';'LA';'ME';'MD';'MA';'MI';'MN';'MS';'MO';...
              'MT';'NE';'NV';'NH';'NJ';'NM';'NY';'NC';'ND';'OH';'OK';'OR';'PA';'RI';'SC';'SD';'TN';'TX';'UT';'VT';'VA';'WA';'WV';'WI';'WY';};        
Fleet.BuildList(nFleetBuild+1,1) = cellstr(get(handles.editBuildName,'string'));
Fleet.Size(nFleetBuild+1,1) = str2double(get(handles.editAvgLoad,'String'));
Fleet.State(nFleetBuild+1,1) = stateAbrev(State);
Fleet.BuildType(nFleetBuild+1,1) = buildTypeName(type);
Fleet.Climate(nFleetBuild+1,1) = cellstr(ClimateZone(char(stateAbrev(State)),str2double(get(handles.editZipCode,'String'))));
Fleet.ZipCode(nFleetBuild+1,1) = str2double(get(handles.editZipCode,'String'));
Fleet.SolarSize(nFleetBuild+1,1) = str2double(get(handles.editSolar,'String'));
set(handles.lbFleet,'string',Fleet.BuildList,'value',1)
lbFleet_Callback(handles.lbFleet, eventdata, handles)

% --- Executes on button press in pushbuttonRemove.
function pushbuttonRemove_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonRemove (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Fleet
val = get(handles.lbFleet,'value');
nFleetBuild = length(Fleet.BuildList);
FleetTemp = Fleet;
clear Fleet
Fleet.Name = FleetTemp.Name;
for i = 1:1:nFleetBuild-1
    if i >= val
        j = i+1;
    else j = i;
    end
    Fleet.BuildList(i,1) = FleetTemp.BuildList(j,1);
    Fleet.Size(i,1) = FleetTemp.Size(j,1);
    Fleet.State(i,1) = FleetTemp.State(j,1);
    Fleet.BuildType(i,1) = FleetTemp.BuildType(j,1);
    Fleet.Climate(i,1) = FleetTemp.Climate(j,1);
    
end

% --- Executes on button press in pushbuttonClear.
function pushbuttonClear_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonClear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Fleet
Fleet.BuildList = {};
Fleet.Size = [];
Fleet.State = {};
Fleet.BuildType = {};
Fleet.Climate = {};
Fleet.ZipCode = [];
Fleet.SolarSize = [];
set(handles.lbFleet,'string',Fleet.BuildList,'value',0)

% --- Executes on selection change in popupmenuBuildType.
function popupmenuBuildType_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuBuildType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
val = get(hObject,'Value');
buildSizes = [43 25 130 400 775 90 9 36 260 82 900 173 49 42 225 47];
set(handles.editAvgLoad,'string',num2str(buildSizes(val)));

% --- Executes during object creation, after setting all properties.
function popupmenuBuildType_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuBuildType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in popupmenuState.
function popupmenuState_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuState (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object creation, after setting all properties.
function popupmenuState_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuState (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editZipCode_Callback(hObject, eventdata, handles)
% hObject    handle to editZipCode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object creation, after setting all properties.
function editZipCode_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editZipCode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editSolar_Callback(hObject, eventdata, handles)
% hObject    handle to editSolar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object creation, after setting all properties.
function editSolar_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editSolar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editAvgLoad_Callback(hObject, eventdata, handles)
% hObject    handle to editAvgLoad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function editSaveName_Callback(hObject, eventdata, handles)
% hObject    handle to editSaveName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object creation, after setting all properties.
function editSaveName_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function editAvgLoad_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editAvgLoad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editBuildName_Callback(hObject, eventdata, handles)
% hObject    handle to editBuildName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object creation, after setting all properties.
function editBuildName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editBuildName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function popupmenuFleet_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuFleet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pushbuttonRunFleet.
function pushbuttonRunFleet_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonRunFleet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Project Model_dir Fleet
h=waitbar(0,'Running Web Tool');

%% load analysis type (set System, FixSize, Control & heat recovery)
if get(handles.uipanelSysSize,'SelectedObject') == handles.radiobuttonFixSize
    SizeStrategy = 1;
elseif get(handles.uipanelSysSize,'SelectedObject') == handles.radiobuttonOptSize
    SizeStrategy = 2;
elseif get(handles.uipanelSysSize,'SelectedObject') == handles.radiobuttonBestNPV
    SizeStrategy = 3;
elseif get(handles.uipanelSysSize,'SelectedObject') == handles.radiobuttonLowestCO2
    SizeStrategy = 4;
end

if get(handles.uipanelSystem,'SelectedObject') ~= handles.radiobuttonCurrent
    ProjDir = fullfile(Model_dir, 'project'); 
    load(fullfile(ProjDir,char('DefaultSystem.mat')))  
end
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
ControlDir = fullfile(Model_dir, 'System Library','Control'); 
if get(handles.uipanelControl,'SelectedObject') == handles.radiobuttonBaseload
    load(fullfile(ControlDir,'1_BaseLoad.mat'));
else
    load(fullfile(ControlDir,'4_LoadFollow.mat'))
end
Project.Control = component;

if get(handles.uipanelCHP,'SelectedObject') == handles.radiobuttonHeatRecov
    DistHeat = 1;
else DistHeat = 0;
end

%% state by state default climates/utilities
stateAbrev = {'AL';'AK';'AZ';'AR';'CA';'CO';'CT';'DE';'FL';'GA';'HI';'ID';'IL';'IN';'IA';'KS';'KY';'LA';'ME';'MD';'MA';'MI';'MN';'MS';'MO';...
              'MT';'NE';'NV';'NH';'NJ';'NM';'NY';'NC';'ND';'OH';'OK';'OR';'PA';'RI';'SC';'SD';'TN';'TX';'UT';'VT';'VA';'WA';'WV';'WI';'WY';}; 
stateName = {'Alabama';'Alaska';'Arizona';'Arkansas';'California';'Colorado';'Connecticut';'Delaware';'Florida';'Georgia';
             'Hawaii';'Idaho';'Illinois';'Indiana';'Iowa';'Kansas';'Kentucky';'Louisiana';'Maine';'Maryland';
             'Massachusetts';'Michigan';'Minnesota';'Mississippi';'Missouri';'Montana';'Nebraska';'Nevada';'NewHampshire';'NewJersey';
             'NewMexico';'NewYork';'NorthCarolina';'NorthDakota';'Ohio';'Oklahoma';'Oregon';'Pennsylvania';'RhodeIsland';'SouthCarolina';
             'SouthDakota';'Tennessee';'Texas';'Utah';'Vermont';'Virginia';'Washington';'WestVirginia';'Wisconsin';'Wyoming';};

BuildDir = fullfile(Model_dir, 'System Library','Buildings'); 
SolarDir = fullfile(Model_dir, 'System Library','Solar','solarData'); 
load(fullfile(SolarDir,'GlobalHorizontal.mat'));
load(fullfile(SolarDir,'Azimuth.mat'));
load(fullfile(SolarDir,'Zenith.mat'));
solar = Project.Renewable.Solar;
OriginalSystem = Project.System;

nFleetBuild = length(Fleet.BuildList);
for build = 1:1:nFleetBuild
    
    %load new building profile
    BuildState = char(Fleet.State(build));
    state = find(strcmp(BuildState,stateAbrev));
    climate = char(Fleet.Climate(build));
    name = strcat(char(Fleet.BuildType(build)),'_',climate,'_New2010.mat');

    load(fullfile(BuildDir,char(name)))
    building.Name = Fleet.BuildList(build);
    building.DistHeat = DistHeat;
    building.DistCool = DistHeat;
    scale = (Fleet.Size(build)/mean(component.DemandE));
    building.DemandE = component.DemandE*scale;
    building.DemandC = building.DistCool*component.DemandC*scale;
    building.DemandH = building.DistHeat*component.DemandH*scale;
    building.NonDistH =(1-building.DistHeat)*component.DemandH*scale;
    building.CoolingElectricalLoad =  building.DistCool*component.CoolingElectricalLoad*scale;
    building.CHPtemp = Project.Building(1).CHPtemp; 
    Project.Building = building;
    
    %% load new solar
    if Fleet.SolarSize(build)>0
        Project.Renewable.Solar = solar;
        Project.Renewable.Solar.State=stateName(state);
        scale = Fleet.SolarSize(build)/Project.Renewable.Solar.Sizem2;
        Project.Renewable.Solar.Sizem2 = Fleet.SolarSize(build);
        Project.Renewable.Solar.Size = Project.Renewable.Solar.Size*scale;
        Project.Renewable.Solar.Irrad = GlobalHorizontal(:,state);
        Project.Renewable.Solar.SunAz = Azimuth(:,state);
        Project.Renewable.Solar.SunZen = Zenith(:,state);
        Project.Renewable.Solar.Tilt = OptimalZenith(Project.Renewable.Solar.Irrad,Project.Renewable.Solar.SunAz,Project.Renewable.Solar.SunZen);
    elseif isfield(Project.Renewable,'Solar')
        Project = rmfield(Project,'Renewable');
    end
    Project.System = OriginalSystem;
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
    RESULTS = runAnalyses(Project);

    Fleet.BuildName(build) = Fleet.BuildList(build);
    Fleet.BuildSize(build) = Fleet.Size(build);
    Fleet.BuildState(build) = Fleet.State(build);
    Fleet.BuildType(build) = Fleet.BuildType(build);
    Fleet.BuildClimate(build) = Fleet.Climate(build);
    Fleet.BuildZipCode(build) = Fleet.ZipCode(build);
    Fleet.BuildSolar(build) = Fleet.SolarSize(build);
    Fleet.NPVbaseline(1:5,build) = [RESULTS.costOut.NPVbaselineDemCharges RESULTS.costOut.NPVbaselineUseCharges  RESULTS.costOut.NPVbaselineFuelCost    RESULTS.costOut.NPVbaselineOandM RESULTS.costOut.NPVbaselineFinance;];
    Fleet.NPVdispatch(1:5,build) = [RESULTS.costOut.NPVnewDemCharges RESULTS.costOut.NPVnewUseCharges RESULTS.costOut.NPVnewFuelCost  RESULTS.costOut.NPVnewOandM RESULTS.costOut.NPVnewFinance;];

    Fleet.InstalledCHPcost(build) =RESULTS.costOut.InstalledCHPcost;
    Fleet.DispatchGridkWhTable(build) = sum(RESULTS.Dispatch.Elec);
    Fleet.BaselineGridkWh(build) = sum(RESULTS.Baseline.Elec);
    Fleet.BaselineHeatkWh(build) = sum(RESULTS.Baseline.Heat);
    Fleet.DispatchHeatkWh(build) = sum(RESULTS.eOut.BoilerHeatHour);
    
    [BaselineEmission DispatchEmission] = EmissionsCalculated(RESULTS,Fleet.State(build),1,1);
    
    Fleet.Summary.CO2baseline(build,1:3) = BaselineEmission(1,1:3);
    Fleet.Summary.NOxbaseline(build,1:3) = BaselineEmission(2,1:3);
    Fleet.Summary.SO2baseline(build,1:3) = BaselineEmission(3,1:3);
    Fleet.Summary.CO2dispatch(build,1:3) = DispatchEmission(1,1:3);
    Fleet.Summary.NOxdispatch(build,1:3) = DispatchEmission(2,1:3);
    Fleet.Summary.SO2dispatch(build,1:3) = DispatchEmission(3,1:3);

    waitbar(build/nFleetBuild,h,strcat('Running Analysis',num2str(build),' of ',num2str(nFleetBuild)));

end
close(h)

ResultDir=fullfile(Model_dir, 'results');
name = strcat('Fleet',get(handles.editSaveName,'String'));
filename=fullfile(ResultDir,name);
save(filename,'Fleet')
close(gcf)
%% Visualization                
FleetResult()

function Tilt = OptimalZenith(Irrad,Az,Zen)
Ao =180;
minZen = min(Zen+99*(Zen==0));
meanZen = sum(Zen)/sum(Zen>0);

% Test multiple tilt angles
To = linspace(minZen,meanZen,20);
Col = linspace(1,1,length(To));

Po = (Irrad*Col).*cosd(Zen*Col-linspace(1,1,length(Zen))'*To).*(cosd(Az-Ao)*Col);
AnPow = sum(Po,1);
[MaxPow, I] = max(AnPow);
if I ==1
    I =2;
end
To = linspace(To(I-1),To(I+1),20);
Po = (Irrad*Col).*cosd(Zen*Col-linspace(1,1,length(Zen))'*To).*(cosd(Az-Ao)*Col);
AnPow = sum(Po,1);
[MaxPow, I] = max(AnPow);
Tilt = To(I);

