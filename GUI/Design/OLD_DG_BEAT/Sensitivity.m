function varargout = Sensitivity(varargin)
% SENSITIVITY M-file for Sensitivity.fig
%      SENSITIVITY, by itself, creates a new SENSITIVITY or raises the existing
%      singleton*.
%
%      H = SENSITIVITY returns the handle to a new SENSITIVITY or the handle to
%      the existing singleton*.
%
%      SENSITIVITY('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SENSITIVITY.M with the given input arguments.
%
%      SENSITIVITY('Property','Value',...) creates a new SENSITIVITY or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Sensitivity_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Sensitivity_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Sensitivity

% Last Modified by GUIDE v2.5 25-May-2014 12:22:09

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Sensitivity_OpeningFcn, ...
                   'gui_OutputFcn',  @Sensitivity_OutputFcn, ...
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


% --- Executes just before Sensitivity is made visible.
function Sensitivity_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Sensitivity (see VARARGIN)

% Choose default command line output for Sensitivity
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Sensitivity wait for user response (see UIRESUME)
% uiwait(handles.figure1);
global  Project Model_dir

BuildType = {'Restaurant: full-service (sit down)';'Restaurant: quick-service (fast food)';'School: primary school';'School: secondary school';'Office: large office';'Office: medium office';'Office: small office';'Mid-rise apartment building';'Hospitality: large hotel';'Hospitality: small hotel/motel';'Health care: large hospital';'Health care: outpatient facility';'Retail: big-box, standalone retail store';'Retail: retail store located in a strip mall';'Retail: supermarket';'Unrefrigerated warehouse';'Multiple'};
Climate = {'Miami (ASHRAE 1A)';'Houston (ASHRAE 2A)';'Phoenix (ASHRAE 2B)';'Atlanta (ASHRAE 3A)';'Las Vegas (ASHRAE 3B-Inland)';'Los Angeles (ASHRAE 3B-Coast)';'San Francisco (ASHRAE 3C)';'Baltimore (ASHRAE 4A)';'Albuquerque (ASHRAE 4B)';'Seattle (ASHRAE 4C)';'Chicago (ASHRAE 5A)';'Boulder (ASHRAE 5B)';'Minneapolis (ASHRAE 6A)';'Helena, MT (ASHRAE 6B)';'Duluth, MN (ASHRAE 7)';'Fairbanks, AK (ASHRAE 8)';};
Vintage = {'2010 construction (ASHRAE 90.1-2010)';'2007 construction (ASHRAE 90.1-2007)';'2004 construction 90.1-2004';'“Post-1980” construction (ASHRAE 90.1-1989)';'“Pre-1980” construction';};

Control_dir=fullfile(Model_dir, 'System Library', 'Control');
ControlFiles=dir(fullfile(Control_dir,'*.mat'));
ControlList=strrep({ControlFiles.name},'.mat','');
name = char(Project.Control.Name);
ControlI = find(strcmp(name,ControlList));
Grid_dir=fullfile(Model_dir, 'System Library', 'Grid');
GridFiles=dir(fullfile(Grid_dir,'*.mat'));
GridList=strrep({GridFiles.name},'.mat','');
name = char(Project.Utilities.Grid.Name);
GridI = find(strcmp(name,GridList));

set(handles.popupmenuControl,'string',ControlList,'value',ControlI(1))
set(handles.popupmenuUtility,'string',GridList,'value',GridI(1))

vintageName = {'New2010';'New2007';'New2004';'Post1980';'Pre1980';};
climateName = {'_1A_'; '_2A_'; '_2B_'; '_3A_'; '_3B_'; '_3B-Coast_'; '_3C_'; '_4A_'; '_4B_'; '_4C_'; '_5A_'; '_5B_'; '_6A_'; '_6B_'; '_7_'; '_8_';};
buildTypeName = {'SDRest'; 'FFRest'; 'Sch-pri'; 'Sch-sec'; 'LgOff'; 'MdOff'; 'SmOff'; 'MRapt'; 'LgHotel'; 'SmHotel'; 'Hospital'; 'OutP'; 'Retail'; 'StMall'; 'SMarket'; 'ware';};

if strcmp(Project.Building.Name,'Multiple') == 1
    name = char(Project.Building.NameList(1));
else name = Project.Building.Name;
end

if strcmp(class(name), 'cell')
    name=name{1};
end

ind =strfind(name,'_');
BuildVal(1) = find(strcmp(name(1:ind(1)-1),buildTypeName));
BuildVal(2) = find(strcmp(name(ind(1):ind(2)),climateName));
BuildVal(3) = find(strcmp(name(ind(2)+1:end),vintageName));

set(handles.popupmenuBuilding,'string',BuildType,'value',BuildVal(1))
set(handles.popupmenuClimate,'string',Climate,'value',BuildVal(2))
set(handles.popupmenuVintage,'string',Vintage,'value',BuildVal(3))
set(handles.checkboxDistHeat,'value',Project.Building.DistHeat(1))
set(handles.checkboxDistCool,'value',Project.Building.DistCool(1))

set(handles.uipanelSysSize,'SelectedObject', handles.radiobuttonOptSize)
set(handles.editFileName,'String',Project.Name);

set(handles.uipanelEcon,'SelectedObject', handles.Stack_Cost)
set(handles.editMin,'string',Project.Economic.InstallCost)
set(handles.editMax,'string',Project.Economic.InstallCost)
set(handles.textNomValue,'string',Project.Economic.InstallCost)
set(handles.editSteps,'string',10)


fields=fieldnames(Project.System);
str={};
for i=1:length(fields)
    if isstruct(Project.System.(fields{i}))
          str{end+1}=fields{i};
    end
end
if max(strcmp(str,'Chiller'))>0
    str2 = Project.System.Chiller(:).ChillType;
    if max(strcmp(str2,'electric'))==0
        set(handles.ElecChill,'enable','off')
    end
    if max(strcmp(str2,'absorption'))==0
        set(handles.AbsorbChill,'enable','off')
    end
end
if max(strcmp(str,'TES'))>0
    str2 = Project.System.TES(:).TEStype;
    if max(strcmp(str2,'hot'))==0
        set(handles.Hot_Storage,'enable','off')
    end
    if max(strcmp(str2,'cold'))==0
        set(handles.Cold_Storage,'enable','off')
    end
end
if max(strcmp(str,'Battery'))==0
    set(handles.Battery_Size,'enable','off')
end
if max(strcmp(str,'PV'))==0
    set(handles.SolarPV,'enable','off')
end


% --- Outputs from this function are returned to the command line.
function varargout = Sensitivity_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on selection change in popupmenuBuilding.
function popupmenuBuilding_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuBuilding (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function popupmenuBuilding_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuBuilding (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function popupmenuClimate_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuClimate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function popupmenuClimate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuClimate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function popupmenuVintage_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuVintage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function popupmenuVintage_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuVintage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function popupmenuControl_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuControl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function popupmenuControl_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuControl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function popupmenuUtility_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuUtility (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function popupmenuUtility_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuUtility (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function checkboxDistCool_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxDistCool (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function checkboxDistHeat_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxDistHeat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function checkboxAllBuild_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxAllBuild (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(handles.checkboxAllBuild,'value')
    set(handles.checkboxAllControl,'value',0)
    set(handles.checkboxAllVintage,'value',0)
    set(handles.checkboxAllClimate,'value',0)
    set(handles.checkboxAllUtility,'value',0)
end

function checkboxAllClimate_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxAllClimate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(handles.checkboxAllClimate,'value')
    set(handles.checkboxAllControl,'value',0)
    set(handles.checkboxAllVintage,'value',0)
    set(handles.checkboxAllUtility,'value',0)
    set(handles.checkboxAllBuild,'value',0)
end

function checkboxAllVintage_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxAllVintage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(handles.checkboxAllVintage,'value')
    set(handles.checkboxAllControl,'value',0)
    set(handles.checkboxAllUtility,'value',0)
    set(handles.checkboxAllClimate,'value',0)
    set(handles.checkboxAllBuild,'value',0)
end

function checkboxAllControl_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxAllControl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(handles.checkboxAllControl,'value')
    set(handles.checkboxAllUtility,'value',0)
    set(handles.checkboxAllVintage,'value',0)
    set(handles.checkboxAllClimate,'value',0)
    set(handles.checkboxAllBuild,'value',0)
end

function checkboxAllUtility_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxAllUtility (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(handles.checkboxAllUtility,'value')
    set(handles.checkboxAllControl,'value',0)
    set(handles.checkboxAllVintage,'value',0)
    set(handles.checkboxAllClimate,'value',0)
    set(handles.checkboxAllBuild,'value',0)
end

% --- Executes when selected object is changed in uipanelSysSize.
function uipanelSysSize_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanelSysSize 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object

function editMin_Callback(hObject, eventdata, handles)
% hObject    handle to editMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object creation, after setting all properties.
function editMin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function editMax_Callback(hObject, eventdata, handles)
% hObject    handle to editMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object creation, after setting all properties.
function editMax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function editSteps_Callback(hObject, eventdata, handles)
% hObject    handle to editSteps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object creation, after setting all properties.
function editSteps_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editSteps (see GCBO)
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

% --- Executes when selected object is changed in uipanelStudy.
function uipanelStudy_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanelStudy 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
if get(handles.uipanelStudy,'SelectedObject') == handles.radiobuttonEcon
    uipanelEcon_SelectionChangeFcn(handles.uipanelEcon, eventdata, handles)
elseif get(handles.uipanelStudy,'SelectedObject') == handles.radiobuttonSys
    uipanelSystem_SelectionChangeFcn(handles.uipanelSystem, eventdata, handles)
end

% --- Executes when selected object is changed in uipanelChillType.
function uipanelEcon_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanelChillType 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
global Project
if get(handles.uipanelStudy,'SelectedObject') == handles.radiobuttonEcon
    sel = get(handles.uipanelEcon,'SelectedObject');
switch get(sel,'Tag')
    case 'Stack_Cost'
        set(handles.editMin,'string',Project.Economic.InstallCost)
        set(handles.editMax,'string',Project.Economic.InstallCost)
        set(handles.textNomValue,'string',Project.Economic.InstallCost)
    case 'Incentive'
        set(handles.editMin,'string',Project.Economic.Incentive)
        set(handles.editMax,'string',Project.Economic.Incentive)
        set(handles.textNomValue,'string',Project.Economic.Incentive)
    case 'SellBack'
        set(handles.editMin,'string',Project.Utilities.Grid.SellBackRate)
        set(handles.editMax,'string',Project.Utilities.Grid.SellBackRate)
        set(handles.textNomValue,'string',Project.Utilities.Grid.SellBackRate)
    case 'Lifespan'
        set(handles.editMin,'string',Project.Economic.LifeYrs)
        set(handles.editMax,'string',Project.Economic.LifeYrs)
        set(handles.textNomValue,'string',Project.Economic.LifeYrs)
    case 'Inflation'
        set(handles.editMin,'string',Project.Economic.Inflation)
        set(handles.editMax,'string',Project.Economic.Inflation)
        set(handles.textNomValue,'string',Project.Economic.Inflation)
    case 'Interest'
        set(handles.editMin,'string',Project.Economic.Interest)
        set(handles.editMax,'string',Project.Economic.Interest)
        set(handles.textNomValue,'string',Project.Economic.Interest)
    case 'Chiller_Cost'
        set(handles.editMin,'string',Project.Economic.ElecChill)
        set(handles.editMax,'string',Project.Economic.ElecChill)
        set(handles.textNomValue,'string',Project.Economic.ElecChill)
    case 'TES_cost'
        set(handles.editMin,'string',Project.Economic.ColdStore)
        set(handles.editMax,'string',Project.Economic.ColdStore)
        set(handles.textNomValue,'string',Project.Economic.ColdStore)
end
set(handles.editSteps,'string',10) 
end

% --- Executes when selected object is changed in uipanelSystem.
function uipanelSystem_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanelSystem 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
global Project
if get(handles.uipanelStudy,'SelectedObject') == handles.radiobuttonSys  
    sel = get(handles.uipanelSystem,'SelectedObject');
switch get(sel,'Tag')
    case 'CHP_Size'
        set(handles.editMin,'string',Project.System.CHP.SysSize(1,1))
        set(handles.editMax,'string',Project.System.CHP.SysSize(1,1))
        set(handles.textNomValue,'string',Project.System.CHP.SysSize(1,1))
    case 'ElecChill'
        ChillSize = 0;
        for i = 1:1:length(Project.System.Chiller)
            if strcmp(Project.System.Chiller(i).ChillType,'electric')
                ChillSize = ChillSize+Project.System.Chiller(i).SysSize;
            end
        end
        set(handles.editMin,'string',ChillSize)
        set(handles.editMax,'string',ChillSize)
        set(handles.textNomValue,'string',ChillSize)
    case 'AbsorbChill'
        ChillSize = 0;
        for i = 1:1:length(Project.System.Chiller)
            if strcmp(Project.System.Chiller(i).ChillType,'absorption')
                ChillSize = ChillSize+Project.System.Chiller(i).SysSize;
            end
        end
        set(handles.editMin,'string',ChillSize)
        set(handles.editMax,'string',ChillSize)
        set(handles.textNomValue,'string',ChillSize)
    case 'Cold_Storage'
        StorageSize = 0;
        for i = 1:1:length(Project.System.TES)
            if strcmp(Project.System.TES(i).TEStype,'cold')
                StorageSize = StorageSize+Project.System.TES(i).SysSize;
            end
        end
        set(handles.editMin,'string',StorageSize)
        set(handles.editMax,'string',StorageSize)
        set(handles.textNomValue,'string',StorageSize)
    case 'Hot_Storage'
        StorageSize = 0;
        for i = 1:1:length(Project.System.TES)
            if strcmp(Project.System.TES(i).TEStype,'hot')
                StorageSize = StorageSize+Project.System.TES(i).SysSize;
            end
        end
        set(handles.editMin,'string',StorageSize)
        set(handles.editMax,'string',StorageSize)
        set(handles.textNomValue,'string',StorageSize)
    case 'Battery_Size'
        StorageSize = sum(Project.System.Battery(:).SysSize);
        set(handles.editMin,'string',StorageSize)
        set(handles.editMax,'string',StorageSize)
        set(handles.textNomValue,'string',StorageSize)
    case 'SolarPV'
        PVSize = sum(Project.System.PV(:).SysSize);
        set(handles.editMin,'string',PVSize)
        set(handles.editMax,'string',PVSize)
        set(handles.textNomValue,'string',PVSize)
end
set(handles.editSteps,'string',10) 
end


% --- Executes on button press in pushbuttonRun.
function pushbuttonRun_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonRun (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Project Model_dir
ProjectOriginal = Project;
%% load analysis type (set System, FixSize, Control & heat recovery)
if get(handles.uipanelStudy,'SelectedObject')==handles.radiobuttonEcon || get(handles.uipanelStudy,'SelectedObject')==handles.radiobuttonSys
    SizeStrategy = 1; %dont resize
elseif get(handles.uipanelSysSize,'SelectedObject') == handles.radiobuttonFixSize
    SizeStrategy = 1; %dont resize
elseif get(handles.uipanelSysSize,'SelectedObject') == handles.radiobuttonOptSize
    SizeStrategy = 2;
elseif get(handles.uipanelSysSize,'SelectedObject') == handles.radiobuttonBestNPV
    SizeStrategy = 3;
elseif get(handles.uipanelSysSize,'SelectedObject') == handles.radiobuttonLowestCO2
    SizeStrategy = 4;
end
count = 0;
h=waitbar(0,'Running Analysis');
if get(handles.uipanelStudy,'SelectedObject')==handles.radiobuttonApp
    sensitivity.StudyType =1;
    GridList = get(handles.popupmenuUtility,'string');
    ControlList = get(handles.popupmenuControl,'string');
    if get(handles.checkboxAllBuild,'value')
        NumRun = 16;
        sensitivity.xlabel ='Building Type';
        sensitivity.xTicks =get(handles.popupmenuBuilding,'string');
    elseif get(handles.checkboxAllClimate,'value')
        NumRun = 16;
        sensitivity.xlabel ='Climate Type';
        sensitivity.xTicks =get(handles.popupmenuClimate,'string');
    elseif get(handles.checkboxAllVintage,'value')
        NumRun = 5;
        sensitivity.xlabel ='Vintage Type';
        sensitivity.xTicks =get(handles.popupmenuVintage,'string');
    elseif get(handles.checkboxAllControl,'value')
        NumRun = length(ControlList);
        sensitivity.xlabel ='Control Type';
        sensitivity.xTicks =get(handles.popupmenuControl,'string');
    elseif get(handles.checkboxAllUtility,'value')
        NumRun = length(GridList);
        sensitivity.xlabel ='Utility Name';
        sensitivity.xTicks =get(handles.popupmenuUtility,'string');
    end
elseif get(handles.uipanelStudy,'SelectedObject')==handles.radiobuttonEcon || get(handles.uipanelStudy,'SelectedObject')==handles.radiobuttonSys
    if get(handles.uipanelStudy,'SelectedObject')==handles.radiobuttonSys
        sensitivity.StudyType =2;
        sensitivity.xlabel =get(get(handles.uipanelSystem,'SelectedObject'),'Tag');
    else sensitivity.StudyType =3;
        sensitivity.xlabel =get(get(handles.uipanelEcon,'SelectedObject'),'Tag');
    end
    NumRun = ceil(str2double(get(handles.editSteps,'string')));
    Min = floor(str2double(get(handles.editMin,'string')));
    Max = ceil(str2double(get(handles.editMax,'string')));
    Steps = ceil(str2double(get(handles.editSteps,'string')));
    sensitivity.X = linspace(Min,Max,Steps);
    sensitivity.xTicks =cellstr(num2str(linspace(Min,Max,Steps)));
end
loadBuilding(1,handles)
loadControl(1,handles)
loadUtility(1,handles)
for runNum= 1:1:NumRun
    count = count+1;
    if sensitivity.StudyType ==1
        loadBuilding(runNum,handles)
        loadControl(runNum,handles)
        loadUtility(runNum,handles)
    elseif sensitivity.StudyType ==2
        changeSysSize(Min+(runNum-1)*(Max-Min)/(Steps-1),handles)
    elseif sensitivity.StudyType ==3
        loadEcon(runNum,handles)
    end
    
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
    if sensitivity.StudyType ==1 || sensitivity.StudyType ==2 || (sensitivity.StudyType ==3 && runNum ==1)
        RESULT = runAnalyses(Project);
    else RESULT.costOut=FinancialCalcs2(RESULT.Baseline,RESULT.Dispatch,Project.Utilities.Grid,Project.Economic,Project.State,'Commercial',Project.Utilities.NatGas);
    end
    % System Summary
    sensitivity.CHPSize(count) = RESULT.Dispatch.SysSize;
    sensitivity.ChillerSize(count) = RESULT.Dispatch.ChillerSize(1);
    sensitivity.TESsize(count)= RESULT.Dispatch.TESsize;
    sensitivity.BatSize(count)= RESULT.Dispatch.BatterySize;
    sensitivity.AvgDemand(count) = mean(Project.Building.DemandE);

    %Financial Results
    Cost = RESULT.costOut;
    sensitivity.BuildName(count) = Project.Building.NameList(1);
    sensitivity.NPVdispatch(count,1:5) = [Cost.NPVnewDemCharges Cost.NPVnewUseCharges Cost.NPVnewFuelCost Cost.NPVnewOandM Cost.NPVnewFinance];
    sensitivity.NPVbaseline(count,1:5) = [Cost.NPVbaselineDemCharges Cost.NPVbaselineUseCharges Cost.NPVbaselineFuelCost Cost.NPVbaselineOandM Cost.NPVbaselineFinance];
    sensitivity.Payback(count) = Cost.Payback;
    sensitivity.Year1baseline(count,1:3) = Cost.Year1baseCharges;
    sensitivity.Year1dispatch(count,1:3) = Cost.Year1dispatchCharges;
    sensitivity.SelfGen(count) = RESULT.eOut.SelfGen;
    sensitivity.CapacitykW(count) = sensitivity.CHPSize(count);
    sensitivity.GenkWh(count) = RESULT.eOut.Total_Electricity_produced_kWh;
    sensitivity.CapacityFactor(count) = sensitivity.GenkWh(count)/(8760*sensitivity.CapacitykW(count));

      
waitbar((count/NumRun),h,strcat('Running Analysis  ',num2str(count),'  of  ',num2str(NumRun)));
end
close(h)

FCModel_dir=fullfile(Model_dir, 'results');  %strrep(which('NREL_FCModel.m'),'\main\NREL_FCModel.m','\results');
name = strcat('sensitivity',get(handles.editFileName,'String'));
filename=fullfile(FCModel_dir,name);
save(filename,'sensitivity')
close(gcf)
%% Visualization
SensitivityResult()
Project = ProjectOriginal;
clear ProjectOriginal


function loadBuilding(runNum,handles)
global Project Model_dir
Project.Building.DistHeat = get(handles.checkboxDistHeat,'value');
Project.Building.DistCool = get(handles.checkboxDistCool,'value');
type=get(handles.popupmenuBuilding, 'value');
climate=get(handles.popupmenuClimate, 'value');
vintage=get(handles.popupmenuVintage, 'value');
if get(handles.checkboxAllBuild,'value')
    type = runNum;
elseif get(handles.checkboxAllClimate,'value')
    climate = runNum;
elseif get(handles.checkboxAllVintage,'value')
    vintage = runNum;
end

vintageName = {'New2010';'New2007';'New2004';'Post1980';'Pre1980';};
climateName = {'_1A_'; '_2A_'; '_2B_'; '_3A_'; '_3B_'; '_3B-Coast_'; '_3C_'; '_4A_'; '_4B_'; '_4C_'; '_5A_'; '_5B_'; '_6A_'; '_6B_'; '_7_'; '_8_';};
buildTypeName = {'SDRest'; 'FFRest'; 'Sch-pri'; 'Sch-sec'; 'LgOff'; 'MdOff'; 'SmOff'; 'MRapt'; 'LgHotel'; 'SmHotel'; 'Hospital'; 'OutP'; 'Retail'; 'StMall'; 'SMarket'; 'ware';};
buildDir = fullfile(Model_dir, 'System Library', 'Buildings');  %strrep(which('NREL_FCModel.m'),'\main\NREL_FCModel.m','\component library\Building');
name = char(strcat(buildTypeName(type),climateName(climate),vintageName(vintage)));
load(strcat(buildDir,filesep,name))
Project.Building.DistHeat = get(handles.checkboxDistHeat,'value');
Project.Building.DistCool = get(handles.checkboxDistCool,'value');
Project.Building.DemandE = component.DemandE;
Project.Building.DemandC = Project.Building.DistCool*component.DemandC;
Project.Building.DemandH = Project.Building.DistHeat*component.DemandH;
Project.Building.NonDistH = (1-Project.Building.DistHeat)*component.DemandH;
Project.Building.CoolingElectricalLoad = Project.Building.DistCool*component.CoolingElectricalLoad;
Project.Building.NameList(1) = cellstr(name);  


function loadUtility(runNum,handles)
global Project Model_dir
GridList = get(handles.popupmenuUtility,'string');
val = get(handles.popupmenuUtility,'value');
if get(handles.checkboxAllUtility,'value')
    name = char(GridList(runNum));
else name = char(GridList(val));
end
Util_dir = fullfile(Model_dir,'System Library','Grid'); %strrep(which('NREL_FCModel.m'),'\main\NREL_FCModel.m','\component library\Grid');
load(strcat(Util_dir,filesep,name))
Project.Utilities.Grid = component;

function loadControl(runNum,handles)
global Project Model_dir
ControlList = get(handles.popupmenuControl,'string');
val = get(handles.popupmenuControl,'value');
if get(handles.checkboxAllControl,'value')
    name = char(ControlList(runNum));
else name = char(ControlList(val));
end
Control_dir = fullfile(Model_dir,'System Library' , 'Control'); % strrep(which('NREL_FCModel.m'),'\main\NREL_FCModel.m','\component library\Control');
load(strcat(Control_dir,filesep,name))
Project.Control = component;

function loadEcon(runNum,handles)
global Project
Min = floor(str2double(get(handles.editMin,'string')));
Max = ceil(str2double(get(handles.editMax,'string')));
Steps = ceil(str2double(get(handles.editSteps,'string')));
if get(handles.uipanelEcon,'SelectedObject')==handles.Stack_Cost
    Project.Economic.InstallCost = Min+(runNum-1)*(Max-Min)/(Steps-1);
elseif get(handles.uipanelEcon,'SelectedObject')==handles.Incentive
    Project.Economic.Incentive = Min+(runNum-1)*(Max-Min)/(Steps-1);
elseif get(handles.uipanelEcon,'SelectedObject')==handles.SellBack
    Project.Utilities.Grid.SellBackRate = Min+(runNum-1)*(Max-Min)/(Steps-1); 
elseif get(handles.uipanelEcon,'SelectedObject')==handles.Lifespan
    Project.Economic.LifeYrs = Min+(runNum-1)*(Max-Min)/(Steps-1);
elseif get(handles.uipanelEcon,'SelectedObject')==handles.Inflation
    Project.Economic.Inflation = Min+(runNum-1)*(Max-Min)/(Steps-1);
elseif get(handles.uipanelEcon,'SelectedObject')==handles.Interest
    Project.Economic.Interest = Min+(runNum-1)*(Max-Min)/(Steps-1);
elseif get(handles.uipanelEcon,'SelectedObject')==handles.Chiller_Cost
    Project.Economic.ElecChill = Min+(runNum-1)*(Max-Min)/(Steps-1);
elseif get(handles.uipanelEcon,'SelectedObject')==handles.TES_cost
    Project.Economic.ColdStore = Min+(runNum-1)*(Max-Min)/(Steps-1);
end

function changeSysSize(NewSize,handles)
global Project
sel = get(handles.uipanelSystem,'SelectedObject');
switch get(sel,'Tag')
    case 'CHP_Size'
        Project.System.CHP.SysSize(1,1)=NewSize;
    case 'ElecChill'
        OldSize = 0;
        for i = 1:1:length(Project.System.Chiller)
            if strcmp(Project.System.Chiller(i).ChillType,'electric')
                OldSize = OldSize+Project.System.Chiller(i).SysSize;
            end
        end
        for i = 1:1:length(Project.System.Chiller)
            if strcmp(Project.System.Chiller(i).ChillType,'electric')
                Project.System.Chiller(i).SysSize = Project.System.Chiller(i).SysSize/OldSize*NewSize;
            end
        end
    case 'AbsorbChill'
        OldSize = 0;
        for i = 1:1:length(Project.System.Chiller)
            if strcmp(Project.System.Chiller(i).ChillType,'absorption')
                OldSize = OldSize+Project.System.Chiller(i).SysSize;
            end
        end
        for i = 1:1:length(Project.System.Chiller)
            if strcmp(Project.System.Chiller(i).ChillType,'absorption')
                Project.System.Chiller(i).SysSize = Project.System.Chiller(i).SysSize/OldSize*NewSize;
            end
        end
    case 'Cold_Storage'
        OldSize = 0;
        for i = 1:1:length(Project.System.TES)
            if strcmp(Project.System.TES(i).TEStype,'cold')
                OldSize = OldSize+Project.System.TES(i).SysSize;
            end
        end
        for i = 1:1:length(Project.System.TES)
            if strcmp(Project.System.TES(i).TEStype,'cold')
                Project.System.TES(i).SysSize = Project.System.TES(i).SysSize/OldSize*NewSize;
            end
        end
    case 'Hot_Storage'
        OldSize = 0;
        for i = 1:1:length(Project.System.TES)
            if strcmp(Project.System.TES(i).TEStype,'hot')
                OldSize = OldSize+Project.System.TES(i).SysSize;
            end
        end
        for i = 1:1:length(Project.System.TES)
            if strcmp(Project.System.TES(i).TEStype,'hot')
                Project.System.TES(i).SysSize = Project.System.TES(i).SysSize/OldSize*NewSize;
            end
        end
    case 'Battery_Size'
        OldSize = sum(Project.System.Battery(:).SysSize);
        Project.System.Battery(:).SysSize = Project.System.Battery(:).SysSize/OldSize*NewSize;
    case 'SolarPV'
        OldSize = sum(Project.System.PV(:).SysSize);
        Project.System.PV(:).SysSize = Project.System.PV(:).SysSize/OldSize*NewSize;
end
