function varargout = MainScreen1(varargin)
% MAINSCREEN1 MATLAB code for MainScreen1.fig
%      MAINSCREEN1, by itself, creates a new MAINSCREEN1 or raises the existing
%      singleton*.
%
%      H = MAINSCREEN1 returns the handle to a new MAINSCREEN1 or the handle to
%      the existing singleton*.
%
%      MAINSCREEN1('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MAINSCREEN1.M with the given input arguments.
%
%      MAINSCREEN1('Property','Value',...) creates a new MAINSCREEN1 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MainScreen1_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MainScreen1_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MainScreen1

% Last Modified by GUIDE v2.5 03-Apr-2017 11:10:48

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MainScreen1_OpeningFcn, ...
                   'gui_OutputFcn',  @MainScreen1_OutputFcn, ...
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


% --- Executes just before MainScreen1 is made visible.
function MainScreen1_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MainScreen1 (see VARARGIN)

% Choose default command line output for MainScreen1
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes MainScreen1 wait for user response (see UIRESUME)
% uiwait(handles.figure1);
global Plant testSystems selectedSystem Model_dir SYSINDEX
Plant.optimptions.method = 'Planning';
if isfield(Plant.Generator,'OpMatA')
    Plant.Generator = rmfield(Plant.Generator,'OpMatA');
end
if isfield(Plant.Generator,'OpMatB')
    Plant.Generator = rmfield(Plant.Generator,'OpMatB');
end
testSystems = [];
selectedSystem = 1;
SYSINDEX = 1;
set(gcf,'Name','Energy Planning Tool 2017.0.1')
popupmenuAxes_Callback(hObject, eventdata, handles)
setupTabs(hObject, handles)
setupSystemSpec(hObject, handles,1)
updateSystemRep(hObject, eventdata, handles)

files = dir(fullfile(Model_dir, 'Plant','*.mat'));
list=strrep({files.name},'.mat','');
set(handles.popupmenuProjectMain,'string',list)
set(handles.popupmenuProjectMain,'value',strmatch(Plant.Name,list,'exact'))
EditSystem(handles)

% --- Outputs from this function are returned to the command line.
function varargout = MainScreen1_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in saveProject.
function saveProject_Callback(hObject, eventdata, handles)
global Model_dir Plant
[f,p]=uiputfile(fullfile(Model_dir,'Plant','PlantNew.mat'),'Save Plant As...');
if f==0; return; end
Plant.Name=strrep(f,'.mat','');
save([p,f],'Plant')


% --- Executes on selection change in popupmenuProjectMain.
function popupmenuProjectMain_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuProjectMain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Model_dir
% Load file that was selected from the popupmenu
projList = get(handles.popupmenuProjectMain,'String');
projName = projList{get(handles.popupmenuProjectMain,'Value')};
projFile = fullfile(Model_dir,'Plant',projName);
load(projFile);
popupmenuAxes_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function popupmenuProjectMain_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuProjectMain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% Contains setup for tabs
function setupTabs(hObject, handles)
%% Main Tabs
% Assumptions:
% 1. Tags of main tab static text boxes are of form, 'MainTab1',
% 'MainTab2', etc.
% 2. Tags of main tab panels are of form, 'uipanelMain1', 'uipanelMain2',
% etc.
TabText = {'Main Window';'Building Spec';'System Spec';'Control Spec';};
for i = 1:length(TabText)
    j = num2str(i);
    % panel management
    set(handles.(strcat('MainTab',j)),'Units','normalized','String',TabText{i});
    set(handles.(strcat('uipanelMain',j)),'Units','normalized','BorderType','none')
    if i ==1
        pan1pos = get(handles.uipanelMain1,'Position');
    else
        set(handles.(strcat('uipanelMain',j)),'Position',pan1pos)
        set(handles.(strcat('uipanelMain',j)),'Visible','off')
    end
end
%% Sub-Tabs
% Assumptions:
% 1. Sub-tabs are only implemented on the first main tab.
% 2. Tags of sub-tab static text boxes are of form, 'mainSubTab1',
% 'mainSubTab2', etc.
% 3. Tags of main tab panels are of form, 'uipanelMainSub1',
% 'uipanelMainSub2', etc.
subTabText = {'Cost Spec';'Building Spec';'System Spec';};
for i = 1:length(subTabText)
    j = num2str(i);
    % panel management
    set(handles.(strcat('mainSubTab',j)),'Units','normalized','String',subTabText{i});
    set(handles.(strcat('uipanelMainSub',j)),'Units','normalized','BorderType','none')
    if i ==1
        pan1pos = get(handles.uipanelMainSub1,'Position');
    else
        set(handles.(strcat('uipanelMainSub',j)),'Position',pan1pos)
        set(handles.(strcat('uipanelMainSub',j)),'Visible','off')
    end
end


% Main tabs callback
function mainTab_Callback(hObject, eventdata, handles)
n = get(hObject,'Tag');
n = n(end);
m = [];
i = 1;
% Find out which tab is currently selected
while isempty(m) && isfield(handles,strcat('uipanelMain',num2str(i)))
    if strcmp(get(handles.(strcat('uipanelMain',num2str(i))),'Visible'),'on')
        m = i;
    else
        i = i+1;
    end
end
m = num2str(m);

% CRUCIAL IN NEXT 3 STEPS: m, then n.
% Change color
bColor = get(handles.(strcat('MainTab',m)),'BackgroundColor');
set(handles.(strcat('MainTab',m)),'BackgroundColor',max(0,bColor-.1))
bColor = get(handles.(strcat('MainTab',n)),'BackgroundColor');
set(handles.(strcat('MainTab',n)),'BackgroundColor',min(1,bColor+.1))

% Change dimensions
pos = get(handles.(strcat('MainTab',m)),'Position');
set(handles.(strcat('MainTab',m)),'Position',[pos(1),pos(2),pos(3),pos(4)-.003])
pos = get(handles.(strcat('MainTab',n)),'Position');
set(handles.(strcat('MainTab',n)),'Position',[pos(1),pos(2),pos(3),pos(4)+.003])

% Change visibility
set(handles.(strcat('uipanelMain',m)),'Visible','off')
set(handles.(strcat('uipanelMain',n)),'Visible','on')

% % Set up System Spec component diagram
% if n == '3'
%     updateSystemRep(hObject, eventdata, handles)
% end


%% Tab 1 functions

% --- Executes on button press in System1.
function System1_Callback(hObject, eventdata, handles)
%% go to system spec tab to edit


% --- Executes on button press in System2.
function System2_Callback(hObject, eventdata, handles)
%% go to system spec tab to edit


% --- Executes on button press in System3.
function System3_Callback(hObject, eventdata, handles)
%% go to system spec tab to edit

% --- Executes on button press in checkboxSys1.
function checkboxSys1_Callback(hObject, eventdata, handles)

% --- Executes on button press in checkboxSys2.
function checkboxSys2_Callback(hObject, eventdata, handles)

% --- Executes on button press in checkboxSys3.
function checkboxSys3_Callback(hObject, eventdata, handles)


% --- Executes on button press in pushbuttonEDC.
function pushbuttonEDC_Callback(hObject, eventdata, handles)
close
DISPATCH

% Sub-tabs callback
function subTab_Callback(hObject, eventdata, handles)
n = get(hObject,'Tag');
n = n(end);
m = [];
i = 1;
% Find out which tab is currently selected
while isempty(m) && isfield(handles,strcat('uipanelMainSub',num2str(i)))
    if strcmp(get(handles.(strcat('uipanelMainSub',num2str(i))),'Visible'),'on')
        m = i;
    else i = i+1;
    end
end
m = num2str(m);


% CRUCIAL IN NEXT 3 STEPS: set m, then n.
% Change color
bColor = get(handles.(strcat('mainSubTab',m)),'BackgroundColor');
set(handles.(strcat('mainSubTab',m)),'BackgroundColor',max(0,bColor-.1))
bColor = get(handles.(strcat('mainSubTab',n)),'BackgroundColor');
set(handles.(strcat('mainSubTab',n)),'BackgroundColor',min(1,bColor+.1))

% Change dimensions
pos = get(handles.(strcat('mainSubTab',m)),'Position');
set(handles.(strcat('mainSubTab',m)),'Position',[pos(1),pos(2),pos(3),pos(4)-.003])
pos = get(handles.(strcat('mainSubTab',n)),'Position');
set(handles.(strcat('mainSubTab',n)),'Position',[pos(1),pos(2),pos(3),pos(4)+.003])

% Change visibility
set(handles.(strcat('uipanelMainSub',m)),'Visible','off')
set(handles.(strcat('uipanelMainSub',n)),'Visible','on')


% --- Executes on selection change in popupmenuAxes.
function popupmenuAxes_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuAxes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Plant
% Plot Electric Demand for now...
% Default values 
Ts = 8760/length(Plant.Data.Demand.E);
zoom = 1;
startDate = 1;
month = [0 31 28 31 30 31 30 31 31 30 31 30 31];
monthDays = [0 31 59 90 120 151 181 212 243 273 304 334 365];
monthLabel = ['January  '; 'February '; 'March    '; 'April    '; 'May      '; 'June     '; 'July     '; 'August   '; 'September'; 'October  ';'November ' ;'December ';];
day1 = 1;
lastDay = 365;
plotAxis = datenum([2014*ones(12,1) (1:12)' 1*ones(12,1) 0*ones(12,1) 0*ones(12,1) 0*ones(12,1)]);
axes(handles.axesMain)
xlabel('Date')
X = datenum(2014,1,1,0,0,0)+((day1-1)+(Ts/24):(Ts/24):lastDay);
Y = Plant.Data.Demand.E(1+24/Ts*(day1-1):24/Ts*lastDay);
plot(X,Y)
ylabel('Demand (kW)')
set(gca,'xtick',plotAxis)
datetick('x','mmmdd','keepticks')


% --- Executes during object creation, after setting all properties.
function popupmenuAxes_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% Sub-tab Cost Spec functions
function CommitEdits_Callback(hObject, eventdata, handles)

function editFixedCost_Callback(hObject, eventdata, handles)

function editFixedCost_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editPeakRate_Callback(hObject, eventdata, handles)

function editPeakRate_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editOffPeakRate_Callback(hObject, eventdata, handles)

function editOffPeakRate_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editPeakHours_Callback(hObject, eventdata, handles)

function editPeakHours_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editDemandCharge_Callback(hObject, eventdata, handles)

function editDemandCharge_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editSellBack_Callback(hObject, eventdata, handles)

function editSellBack_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editGasCost_Callback(hObject, eventdata, handles)

function editGasCost_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editOthrFixedCosts_Callback(hObject, eventdata, handles)

function editOthrFixedCosts_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% Sub-tab Building Spec functions (call same functions as Building Spec Main Tab)
function sliderArea2_Callback(hObject, eventdata, handles)

function sliderArea2_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function sliderWindowWall2_Callback(hObject, eventdata, handles)

function sliderWindowWall2_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function sliderLighting2_Callback(hObject, eventdata, handles)

function sliderLighting2_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function sliderPlugLoad2_Callback(hObject, eventdata, handles)

function sliderPlugLoad2_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function sliderDensity2_Callback(hObject, eventdata, handles)

function sliderDensity2_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function sliderUseHours2_Callback(hObject, eventdata, handles)

function sliderUseHours2_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function popupmenuLighting_Callback(hObject, eventdata, handles)

function popupmenuLighting_CreateFcn(hObject, eventdata, handles)

%% Sub-tab System Spec functions
function sliderSize_Callback(hObject, eventdata, handles)

function sliderSize_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function sliderMaintenance_Callback(hObject, eventdata, handles)

function sliderMaintenance_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function popupmenuComponent_Callback(hObject, eventdata, handles)

function popupmenuComponent_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% Tab 2 functions
function popupmenuLocation_Callback(hObject, eventdata, handles)

function popupmenuLocation_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function popupmenuBuildType_Callback(hObject, eventdata, handles)

function popupmenuBuildType_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function popupmenuBuildAge_Callback(hObject, eventdata, handles)

function popupmenuBuildAge_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function sliderArea_Callback(hObject, eventdata, handles)

function sliderArea_CreateFcn(hObject, eventdata, handles)
function sliderWindowWall_Callback(hObject, eventdata, handles)

function sliderWindowWall_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function sliderLighting_Callback(hObject, eventdata, handles)

function sliderLighting_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function sliderPlugLoad_Callback(hObject, eventdata, handles)

function sliderPlugLoad_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function sliderDensity_Callback(hObject, eventdata, handles)

function sliderDensity_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function sliderUseHours_Callback(hObject, eventdata, handles)

function sliderUseHours_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

%% Tab 3 functions
function setupSystemSpec(hObject, handles,sysnum)
global Plant Model_dir
%set up the system spec tab
%give the buttons symbols
buttonIms = {'mGT', 'ICE', 'SOFC', 'MCFC', 'PEM', 'SolarPV', 'SolarThermal', 'Chiller', ...
    'AbChiller','AirHeater','WaterHeater', 'ColdStor', 'HotStor', 'Battery'}; % need images for: 'SolarStirling', 'Wind','HighTempStor
for i = 1:1:length(buttonIms)
    button = buttonIms{i};
    imfile = fullfile(Model_dir,'GUI','Graphics',strcat(button, '.png'));
    [x,map] = imread(imfile);
    I2 = imresize(x,[42 70]);
    set(handles.(button), 'cdata', I2)
end
%%load all current components into the GUI and the plant
% --- Executes on button press in System1Comps.
function System1Comps_Callback(hObject, eventdata, handles)
switchsys(1,handles)

function System2Comps_Callback(hObject, eventdata, handles)
switchsys(2,handles)

function System3Comps_Callback(hObject, eventdata, handles)
switchsys(3,handles)

function switchsys(currentSys,handles)
global selectedSystem
if get(handles.(strcat('System',num2str(currentSys),'Comps')),'Value')
    set(handles.(strcat('System',num2str(selectedSystem),'Comps')),'Value',0)
    offcolor = get(handles.(strcat('System',num2str(currentSys),'Comps')),'BackgroundColor');
    oncolor = get(handles.(strcat('System',num2str(selectedSystem),'Comps')),'BackgroundColor');
    set(handles.(strcat('System',num2str(currentSys),'Comps')),'BackgroundColor',oncolor)
    set(handles.(strcat('System',num2str(selectedSystem),'Comps')),'BackgroundColor',offcolor)
    setupSystemSpec([],handles,currentSys)
    selectedSystem = currentSys;
end

function pushbuttonAbChillerInSys_Callback(hObject, eventdata, handles)
SetupSystem(hObject,handles)

function pushbuttonACDC_Callback(hObject, eventdata, handles)
SetupSystem(hObject,handles)

function pushbuttonBatteryInSys_Callback(hObject, eventdata, handles)
SetupSystem(hObject,handles)

function pushbuttonChillerInSys_Callback(hObject, eventdata, handles)
SetupSystem(hObject,handles)

function pushbuttonCoolingDemands_Callback(hObject, eventdata, handles)
SetupSystem(hObject,handles)

function pushbuttonFuelCell_Callback(hObject, eventdata, handles)
SetupSystem(hObject,handles)

function pushbuttonGrid_Callback(hObject, eventdata, handles)
SetupSystem(hObject,handles)

function pushbuttonHeaterInSys_Callback(hObject, eventdata, handles)
SetupSystem(hObject,handles)

function pushbuttonHeatingDemands_Callback(hObject, eventdata, handles)
SetupSystem(hObject,handles)

function pushbuttonHotWaterDemands_Callback(hObject, eventdata, handles)
SetupSystem(hObject,handles)

function pushbuttonICE_mGT_Callback(hObject, eventdata, handles)
SetupSystem(hObject,handles)

function pushbuttonSolarPVInSys_Callback(hObject, eventdata, handles)
SetupSystem(hObject,handles)

function pushbuttonSolarSterlingInSys_Callback(hObject, eventdata, handles)
SetupSystem(hObject,handles)

function pushbuttonSolarThermalInSys_Callback(hObject, eventdata, handles)
SetupSystem(hObject,handles)

function pushbuttonTES1_Callback(hObject, eventdata, handles)
SetupSystem(hObject,handles)

function pushbuttonTES2_Callback(hObject, eventdata, handles)
SetupSystem(hObject,handles)

function pushbuttonTES3_Callback(hObject, eventdata, handles)
SetupSystem(hObject,handles)

function pushbuttonWaterHeaterInSys_Callback(hObject, eventdata, handles)
SetupSystem(hObject,handles)

function pushbuttonWindInSys_Callback(hObject, eventdata, handles)
SetupSystem(hObject,handles)

%% The following need to represent different things for each system selected in 3rd tab
function CompName_Callback(hObject, eventdata, handles)
global Plant
plantindex = get(hObject,'Userdata');
Plant.Generator(plantindex).Name = get(hObject,'String');
ncomps = length(handles.uipanelDem.Children);
child = ncomps-3-(plantindex-3);
set(handles.uipanelDem.Children(child),'String',get(hObject,'String'))

function CompName_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function MaxEff_Callback(hObject, eventdata, handles)
global Plant
plantindex = get(handles.CompName,'Userdata');
Plant.Generator(plantindex).VariableStruct.eff = str2double(get(hObject,'String'));

function MaxEff_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function RampRate_Callback(hObject, eventdata, handles)
global Plant
plantindex = get(handles.CompName,'Userdata');
Plant.Generator(plantindex).VariableStruct.Ramp = str2double(get(hObject,'String'));

function RampRate_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function MinOutput_Callback(hObject, eventdata, handles)
global Plant
plantindex = get(handles.CompName,'Userdata');
Plant.Generator(plantindex).VariableStruct.LB = str2double(get(hObject,'String'));

function MinOutput_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function MaxOutput_Callback(hObject, eventdata, handles)
global Plant
plantindex = get(handles.CompName,'Userdata');
Plant.Generator(plantindex).Size = str2double(get(hObject,'String'));
Plant.Generator(plantindex).VariableStruct.UB = str2double(get(hObject,'String'));

function MaxOutput_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function uitableEffCurve_CellEditCallback(hObject, eventdata, handles)
global Plant SYSINDEX
Outputs = fieldnames(Plant.Generator(SYSINDEX).Output);
nOutput = eventdata.Indices;
newValue = eventdata.NewData;
Plant.Generator(SYSINDEX).Output.(Outputs{nOutput(2)})(nOutput(1)) = newValue;

function saveSystem_Callback(hObject, eventdata, handles)
global selectedSystem Plant testSystems Model_dir
testSystems(selectedSystem).Plant = Plant;
[f,p]=uiputfile(fullfile(Model_dir,'Plant','PlantNew.mat'),'Save Plant As...');
if f==0; return; end
Plant.Name=strrep(f,'.mat','');
save([p,f],'Plant')


%% Library (Tab 3) callbacks
function selectComponent(hObject,selected, handles)
%%Give user option to load system of this type from library
%%and add to or replace in system configuration, or to edit one 
%%then add or replace

% --- Executes on Micro Turbine button press.
function mGT_Callback(hObject, eventdata, handles)
AddSystem(hObject,handles)

function ICE_Callback(hObject, eventdata, handles)
AddSystem(hObject,handles)

function CustomCHP_Callback(hObject, eventdata, handles)
AddSystem(hObject,handles)

function SOFC_Callback(hObject, eventdata, handles)
AddSystem(hObject,handles)

function MCFC_Callback(hObject, eventdata, handles)
AddSystem(hObject,handles)

function PEM_Callback(hObject, eventdata, handles)
AddSystem(hObject,handles)

function CustomFC_Callback(hObject, eventdata, handles)
AddSystem(hObject,handles)

function SolarPV_Callback(hObject, eventdata, handles)
AddSystem(hObject,handles)

function SolarThermal_Callback(hObject, eventdata, handles)
AddSystem(hObject,handles)

function SolarStirling_Callback(hObject, eventdata, handles)
AddSystem(hObject,handles)

function Wind_Callback(hObject, eventdata, handles)
AddSystem(hObject,handles)

function Chiller_Callback(hObject, eventdata, handles)
AddSystem(hObject,handles)

function AbChiller_Callback(hObject, eventdata, handles)
AddSystem(hObject,handles)

function AirHeater_Callback(hObject, eventdata, handles)
AddSystem(hObject,handles)

function WaterHeater_Callback(hObject, eventdata, handles)
AddSystem(hObject,handles)

function ColdStor_Callback(hObject, eventdata, handles)
AddSystem(hObject,handles)

function HighTempStor_Callback(hObject, eventdata, handles)
AddSystem(hObject,handles)

function HotStor_Callback(hObject, eventdata, handles)
AddSystem(hObject,handles)

function ElecStor_Callback(hObject, eventdata, handles)
AddSystem(hObject,handles)

% --- Executes on button press in pushbuttonRemove.
function pushbuttonRemove_Callback(hObject, eventdata, handles)
global Plant SYSINDEX
nG = length(Plant.Generator);
str = strcat(Plant.Generator(SYSINDEX).Type,'.',Plant.Generator(SYSINDEX).Name);
for n = 1:1:length(Plant.Network)
    Plant.Network(n).Equipment = Plant.Network(n).Equipment(~strcmp(str,Plant.Network(n).Equipment));
end
if SYSINDEX ==1
    Plant.Generator = Plant.Generator(2:end);
elseif SYSINDEX == nG
    Plant.Generator = Plant.Generator(1:end-1);
else
    Plant.Generator = [Plant.Generator(1:SYSINDEX-1),Plant.Generator(SYSINDEX+1:end)];
end
SYSINDEX = 1;
updateSystemRep(hObject, eventdata, handles)
EditSystem(handles)

% --- Executes on selection change in Compfuel.
function Compfuel_Callback(hObject, eventdata, handles)
global Plant
plantindex = get(handles.CompName,'Userdata');
fuels = get(hObject,'String');
Plant.Generator(plantindex).Source = fuels{get(hObject,'Value')};

% --- Executes during object creation, after setting all properties.
function Compfuel_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% Tab 4 functions
function popupmenuOptimization_Callback(hObject, eventdata, handles)

function popupmenuOptimization_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function AggressiveOpt_Callback(hObject, eventdata, handles)
if get(handles.AggressiveOpt,'Value')
    set(handles.MedAncilOpt,'Value',0)
    set(handles.RobustAncilOpt,'Value',0)
end

function MedAncilOpt_Callback(hObject, eventdata, handles)
if get(handles.MedAncilOpt,'Value')
    set(handles.AggressiveOpt,'Value',0)
    set(handles.RobustAncilOpt,'Value',0)
end

function RobustAncilOpt_Callback(hObject, eventdata, handles)
if get(handles.RobustAncilOpt,'Value')
    set(handles.AggressiveOpt,'Value',0)
    set(handles.MedAncilOpt,'Value',0)
end

function AutoAncilOpt_Callback(hObject, eventdata, handles)
if get(handles.AutoAncilOpt,'Value')
    set(handles.ManualAncilOpt,'Value',0)
else
    set(handles.ManualAncilOpt,'Value',1)
end

function ManualAncilOpt_Callback(hObject, eventdata, handles)
if get(handles.ManualAncilOpt,'Value')
    set(handles.AutoAncilOpt,'Value',0)
else
    set(handles.AutoAncilOpt,'Value',1)
end


%% Need to automatically creat these control menus for all dispatchable generators
function popupmenuSys1_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function popupmenuSys1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% Need to create extra sliders if there are more than 1 building
function sliderComfort_Callback(hObject, eventdata, handles)

function sliderComfort_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



