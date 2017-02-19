function varargout = NationalSurveyResults(varargin)
% NATIONALSURVEYRESULTS MATLAB code for NationalSurveyResults.fig
%      NATIONALSURVEYRESULTS, by itself, creates a new NATIONALSURVEYRESULTS or raises the existing
%      singleton*.
%
%      H = NATIONALSURVEYRESULTS returns the handle to a new NATIONALSURVEYRESULTS or the handle to
%      the existing singleton*.
%
%      NATIONALSURVEYRESULTS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in NATIONALSURVEYRESULTS.M with the given input arguments.
%
%      NATIONALSURVEYRESULTS('Property','Value',...) creates a new NATIONALSURVEYRESULTS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before NationalSurveyResults_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to NationalSurveyResults_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help NationalSurveyResults

% Last Modified by GUIDE v2.5 26-May-2014 12:23:20

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @NationalSurveyResults_OpeningFcn, ...
                   'gui_OutputFcn',  @NationalSurveyResults_OutputFcn, ...
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


% --- Executes just before NationalSurveyResults is made visible.
function NationalSurveyResults_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to NationalSurveyResults (see VARARGIN)

% Choose default command line output for NationalSurveyResults
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
global Model_dir MapHandle
MapHandle = figure(1);
resultsdir=fullfile(Model_dir, 'results');
x=dir(fullfile(resultsdir,'*.mat'));
listTemp=strrep({x.name},'.mat','');
i = 0;
for j = 1:1:length(listTemp')
    if  strncmp(listTemp(j),'NatSurv',7)
        i = i+1;
        list(i) = strrep(listTemp(j),'NatSurv','');
    end
end
if ~isempty(find(strcmp('Project1',list)));
    currentProject = find(strcmp('Project1',list));
else currentProject=1;
end
set(handles.popupmenuResultsFile,'string',list,'value',currentProject)

%plot options
plots1={'Total Commercial Fleet (GW)'
        'Self Generation Fraction'
        'Average FC parity price ($/kW)'
        'Peak FC parity price ($/kW)'
        'Average Net Present Value Savings (%)'
        'Best Net Present Value Savings (%)'
        'Average GHG emission reduction (%)'};

set(handles.popupmenuMap,'string',plots1,'value',1)
popupmenuResultsFile_Callback(handles.popupmenuResultsFile, eventdata, handles)


% --- Outputs from this function are returned to the command line.
function varargout = NationalSurveyResults_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in popupmenuResultsFile.
function popupmenuResultsFile_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuResultsFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global result Model_dir
str=get(hObject,'string');
val=get(hObject,'value'); 
resultsdir=fullfile(Model_dir, 'results');
load(fullfile(resultsdir, strcat('NatSurv',str{val})))
result = Survey;
result.Econ.MarketPen = 1;

set(handles.editMarketPen,'string',1)
set(handles.sliderMarketPen,'Min',0,'Max',5,'Value',.5)
set(handles.uipanelRegion,'SelectedObject', handles.radiobuttonNational)

%load buildings  lists
BuildType = {'Restaurant: full-service (sit down)';'Restaurant: quick-service (fast food)';'School: primary school';'School: secondary school';'Office: large office';'Office: medium office';'Office: small office';'Mid-rise apartment building';'Hospitality: large hotel';'Hospitality: small hotel/motel';'Health care: large hospital';'Health care: outpatient facility';'Retail: big-box, standalone retail store';'Retail: retail store located in a strip mall';'Retail: supermarket';'Unrefrigerated warehouse';};
buildTypeName = {'SDRest'; 'FFRest'; 'Sch-pri'; 'Sch-sec'; 'LgOff'; 'MdOff'; 'SmOff'; 'MRapt'; 'LgHotel'; 'SmHotel'; 'Hospital'; 'OutP'; 'Retail'; 'StMall'; 'SMarket'; 'ware';};
if strcmp(result.BuildName,'all')==1
    list = BuildType;
    set(handles.pushbuttonAddBuild,'enable','on')
    set(handles.pushbuttonAllBuild,'enable','on')
    set(handles.pushbuttonRemoveBuild,'enable','on')
    set(handles.pushbuttonClearBuild,'enable','on')
elseif length(result.BuildName)==1
    bN = find(strcmp(result.BuildName,buildTypeName));
    if isempty(bN)
        list = result.BuildName;
    else list = BuildType(bN);
    end
    set(handles.pushbuttonAddBuild,'enable','off')
    set(handles.pushbuttonAllBuild,'enable','off')
    set(handles.pushbuttonRemoveBuild,'enable','off')
    set(handles.pushbuttonClearBuild,'enable','off')
else
    for j = 1:1:length(result.BuildName)
        bN = find(strcmp(result.BuildName(j),buildTypeName));
        if isempty(bN)
            list(j) = result.BuildName(j);
        else list(j) = BuildType(bN);
        end
    end
    set(handles.pushbuttonAddBuild,'enable','off')
    set(handles.pushbuttonAllBuild,'enable','off')
    set(handles.pushbuttonRemoveBuild,'enable','off')
    set(handles.pushbuttonClearBuild,'enable','off')
end
set(handles.listboxBuildings, 'string', list, 'Max', length(list), 'Min', 1)
set(handles.listboxPlotBuild, 'string', list, 'Max', length(list), 'Min', 1)

%load state lists
stateName = {'Alabama';'Alaska';'Arizona';'Arkansas';'California';'Colorado';'Connecticut';'Delaware';'Florida';'Georgia';
             'Hawaii';'Idaho';'Illinois';'Indiana';'Iowa';'Kansas';'Kentucky';'Louisiana';'Maine';'Maryland';
             'Massachusetts';'Michigan';'Minnesota';'Mississippi';'Missouri';'Montana';'Nebraska';'Nevada';'New Hampshire';'New Jersey';
             'New Mexico';'New York';'North Carolina';'North Dakota';'Ohio';'Oklahoma';'Oregon';'Pennsylvania';'Rhode Island';'South Carolina';
             'South Dakota';'Tennessee';'Texas';'Utah';'Vermont';'Virginia';'Washington';'West Virginia';'Wisconsin';'Wyoming';};
set(handles.listboxStates, 'string', stateName, 'Max', length(stateName), 'Min', 1)
set(handles.listboxPlotStates, 'string', stateName, 'Max', length(stateName), 'Min', 1)
popupmenuMap_Callback(handles.popupmenuMap, eventdata, handles)


% --- Executes on selection change in popupmenuAxes1.
function popupmenuMap_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuAxes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
val=get(hObject,'value');
str=get(hObject,'string');
currentstr=str{val};
stateAbrev = {'AL';'AK';'AZ';'AR';'CA';'CO';'CT';'DE';'FL';'GA';'HI';'ID';'IL';'IN';'IA';'KS';'KY';'LA';'ME';'MD';'MA';'MI';'MN';'MS';'MO';
              'MT';'NE';'NV';'NH';'NJ';'NM';'NY';'NC';'ND';'OH';'OK';'OR';'PA';'RI';'SC';'SD';'TN';'TX';'UT';'VT';'VA';'WA';'WV';'WI';'WY';}; 
statesAll = get(handles.listboxStates,'string');
states = get(handles.listboxPlotStates,'string');
stateNum = [];
for i = 1:1:length(states)
    stateNum(end+1) = find(strcmp(states(i),statesAll));
end
regionList={'NorthEast';'Midwest';'South';'West';};
region = [3 4 4 3 4 4 1 3 3 3 4 4 2 2 2 2 3 3 1 3 1 2 2 3 2 4 2 4 1 1 4 1 3 2 2 3 4 1 1 3 2 3 3 4 1 3 4 3 2 4];

[A,units] = formatResults(currentstr,handles);
Y = zeros(50,1);
Y(stateNum) = A(stateNum);
% Y(isnan(Y)) = 0;
if length(states)>13
    switch currentstr
        case 'Total Commercial Fleet (GW)'
            Y2 = [Y'.*(region==1); Y'.*(region==2); Y'.*(region==3); Y'.*(region==4);];
        otherwise
            Y2 = [sum(Y'.*(region==1))/sum(region==1) sum(Y'.*(region==2))/sum(region==2) sum(Y'.*(region==3))/sum(region==3) sum(Y'.*(region==4))/sum(region==4)];
    end
    AxisTags =regionList;
    LabelX = ('Census Region');
else 
    Y2 = nonzeros(Y)';
    AxisTags = stateAbrev(stateNum);
    LabelX = ('State');
end
axes(handles.axes1);
hold off
if ndims(Y2) >1
    bar(Y2,'stacked')
else bar(Y2)
end
set(gca,'XTickLabel',AxisTags)
ylabel(units)
xlabel(LabelX)

States2Plot= [];
for i = 1:1:length(states)
    if strcmp(statesAll(2),states(i)) || strcmp(statesAll(11),states(i))
        if length(states)<=2 %omit hawaii and alaska
            States2Plot(end+1) = find(strcmp(states(i),statesAll));
        end
    else
        States2Plot(end+1) = find(strcmp(states(i),statesAll));
    end
end
States2Plot = sort(States2Plot);
states = statesAll(States2Plot);
B = sort(A,'descend');
ColorScaleMax = B(1);
if A(11) == ColorScaleMax
    ColorScaleMax = B(2);
end
ColorScaleMin = min(0,B(end));
if ColorScaleMax<15
    ColorScaleMax = ceil(ColorScaleMax);
    ticks = ceil(ColorScaleMax-ColorScaleMin);
elseif ColorScaleMax<100
    ColorScaleMax = 10*ceil(ColorScaleMax/10);
    ticks = ceil((ColorScaleMax-ColorScaleMin)/100);
else ColorScaleMax = 1000*ceil(ColorScaleMax/1000);
    ticks = ceil((ColorScaleMax-ColorScaleMin)/1000);
end
if ColorScaleMin>-10
    ColorScaleMin = floor(ColorScaleMin);
elseif ColorScaleMin>-100
    ColorScaleMin = 10*floor(ColorScaleMin/10);
else ColorScaleMin = 1000*floor(ColorScaleMin/1000);
end
if ticks<3
    ticks = ceil(ColorScaleMax-ColorScaleMin)/10+1;
end
if    ticks>10
    ticks = 10;
end
SupplyDemand(handles)
NationalMapPlot(Y(States2Plot),states,units,ColorScaleMax,ColorScaleMin,ticks,'parula')
%% color: 'parula', 'hsv','jet', 'spring', 'summer', 'autumn'
%% 'gray'

% --- Executes during object creation, after setting all properties.
function popupmenuMap_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuMap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function [Y,units] = formatResults(currentstr,handles)
global result
load USbuildData

buildNames = get(handles.listboxBuildings, 'string');
stateNames =  get(handles.listboxStates, 'string');
list1 = get(handles.listboxPlotBuild, 'string');
% list2 = get(handles.listboxPlotStates, 'string');
list2 = get(handles.listboxStates, 'string');
builds = zeros(50,length(buildNames));
for i = 1:1:length(list1)
    bN = find(strcmp(list1(i),buildNames));
    for j = 1:1:length(list2)
        sN = find(strcmp(list2(j),stateNames));
        builds(sN,bN) =buildData(sN,bN);
    end
end

TotalBuildLoadkWh = result.GenkWh./result.SelfGen;

TCBF = sum(TotalBuildLoadkWh.*builds,2); %total kWh
SGF = sum(result.GenkWh.*builds,2)./TCBF*100; %demensionless
APP  = sum(result.FCevenCost.*result.FCsize.*builds,2)./sum(result.FCsize.*builds,2); %$/kW
PPP = max(result.FCevenCost,[],2); %$/kW
%NPC Savings
NPCbase= squeeze(sum(result.NPCbaseline,1));
NPCnew = squeeze(sum(result.NPCdispatch,1));
NPVS = (NPCbase-NPCnew)./NPCbase*100; %demensionless
ANPVS = sum(NPVS.*TotalBuildLoadkWh.*builds,2)./TCBF;%demensionless
BNPVS = max(NPVS,[],2); %demensionless
%Emissions
A = size(result.Baseline);
for i = 1:1:A(1)
    for j = 1:1:A(2)
        Base(i,j) = result.Baseline(i,j).CO2;
        DispatchCO2(i,j) = result.Dispatch(i,j).CO2;
    end
end
GHGred = sum((Base-DispatchCO2).*builds,2)./sum(Base.*builds,2)*100;

switch currentstr
    case 'Total Commercial Fleet (GW)'
        units = 'Total Demand (GW)';
        Y = TCBF/8760/1e6;%convert to GW
    case 'Self Generation Fraction'
        units = 'Self Gen % of Building Demand';
        Y = SGF;
    case 'Average FC parity price ($/kW)'
        units = 'FC Parity Price ($/kW)';
        Y =APP;
    case 'Peak FC parity price ($/kW)'
        units = 'FC Parity Price ($/kW)';
        Y = PPP;
    case 'Average Net Present Value Savings (%)'
        units = '% of NPC Saved with DG';
        Y = ANPVS;
    case 'Best Net Present Value Savings (%)'
        units = '% of NPC Saved with DG';
        Y = BNPVS;
    case 'Average GHG emission reduction (%)'
        units = '% GHG emission reduction with DG';
        Y = GHGred;
end

% --- Executes during object creation, after setting all properties.
function uipanelRegion_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uipanelRegion (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% --- Executes when selected object is changed in uipanelRegion.
function uipanelRegion_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanelRegion 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
if get(handles.uipanelRegion,'SelectedObject')== handles.radiobuttonNational
    States2Plot = linspace(1,50,50);
%     States2Plot = [1, linspace(3,10,8), linspace(12,50,39)]; %omit hawaii and alaska
elseif get(handles.uipanelRegion,'SelectedObject')== handles.radiobuttonNE
    States2Plot = [7, 8, 19, 20, 21, 29, 30, 32, 38, 39, 45];
elseif get(handles.uipanelRegion,'SelectedObject')== handles.radiobuttonSouth
    States2Plot = [1, 4, 9, 10, 17, 18, 24, 33, 40, 42, 43, 46, 48];    
elseif get(handles.uipanelRegion,'SelectedObject')== handles.radiobuttonMW
    States2Plot = [13, 14, 15, 16, 22, 23, 25, 27, 34, 35, 36, 41, 49];    
elseif get(handles.uipanelRegion,'SelectedObject')== handles.radiobuttonWC
    States2Plot = [3, 5, 6, 12, 26, 28, 31, 37, 44, 47, 50];  
elseif get(handles.uipanelRegion,'SelectedObject')== handles.radiobuttonAKHI
    States2Plot = [2,11];  
end
states = get(handles.listboxStates,'string');
set(handles.listboxPlotStates, 'string',states(States2Plot));
popupmenuMap_Callback(handles.popupmenuMap, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function popupmenuResultsFile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuResultsFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in listboxBuildings.
function listboxBuildings_Callback(hObject, eventdata, handles)
% hObject    handle to listboxBuildings (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function listboxBuildings_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listboxBuildings (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in listboxPlotBuild.
function listboxPlotBuild_Callback(hObject, eventdata, handles)
% hObject    handle to listboxPlotBuild (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object creation, after setting all properties.
function listboxPlotBuild_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listboxPlotBuild (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pushbuttonAddBuild.
function pushbuttonAddBuild_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonAddBuild (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
list1 = get(handles.listboxBuildings, 'string');
list2 = get(handles.listboxPlotBuild, 'string');
val = get(handles.listboxBuildings, 'Value');
if length(list2) ==1 && strcmp(list2(1),'None')
    list2 = {};
end
for i = 1:1:length(val)
    build = char(list1(val(i)));
    if isempty(list2) || max(strcmp(build,list2))==0
        list2(end+1) = cellstr(build);
    end
end
set(handles.listboxPlotBuild, 'string', list2, 'Max', length(list2), 'Min', 1)
popupmenuMap_Callback(handles.popupmenuMap, eventdata, handles)

% --- Executes on button press in pushbuttonRemoveBuild.
function pushbuttonRemoveBuild_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonRemoveBuild (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
list = get(handles.listboxPlotBuild, 'string');
val = get(handles.listboxPlotBuild, 'Value');
list2={};
for i = 1:1:length(list)
    if max(i==val)==0
        list2(end+1) = list(i);
    end
end
set(handles.listboxPlotBuild, 'string', list2, 'Max', length(list2), 'Min', 1)
popupmenuMap_Callback(handles.popupmenuMap, eventdata, handles)

% --- Executes on button press in pushbuttonAllBuild.
function pushbuttonAllBuild_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonAllBuild (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
list = get(handles.listboxBuildings, 'string');
set(handles.listboxPlotBuild, 'string', list, 'Max', length(list), 'Min', 1)
popupmenuMap_Callback(handles.popupmenuMap, eventdata, handles)

% --- Executes on button press in pushbuttonClearBuild.
function pushbuttonClearBuild_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonClearBuild (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
list = {'None'};
set(handles.listboxPlotBuild, 'string', list, 'Max', length(list), 'Min', 1)

% --- Executes on selection change in listboxStates.
function listboxStates_Callback(hObject, eventdata, handles)
% hObject    handle to listboxStates (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object creation, after setting all properties.
function listboxStates_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listboxStates (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in listboxPlotStates.
function listboxPlotStates_Callback(hObject, eventdata, handles)
% hObject    handle to listboxPlotStates (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object creation, after setting all properties.
function listboxPlotStates_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listboxPlotStates (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pushbuttonAddState.
function pushbuttonAddState_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonAddState (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
list1 = get(handles.listboxStates, 'string');
list2 = get(handles.listboxPlotStates, 'string');
val = get(handles.listboxStates, 'Value');
if length(list2) ==1 && strcmp(list2(1),'None')
    list2 = {};
end
for i = 1:1:length(val)
    state = char(list1(val(i)));
    if isempty(list2) || max(strcmp(state,list2))==0
        list2(end+1) = cellstr(state);
    end
end
set(handles.listboxPlotStates, 'string', list2, 'Max', length(list2), 'Min', 1)
popupmenuMap_Callback(handles.popupmenuMap, eventdata, handles)

% --- Executes on button press in pushbuttonRemoveState.
function pushbuttonRemoveState_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonRemoveState (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
list = get(handles.listboxPlotStates, 'string');
val = get(handles.listboxPlotStates, 'Value');
list2={};
for i = 1:1:length(list)
    if max(i==val)==0
        list2(end+1) = list(i);
    end
end
set(handles.listboxPlotStates, 'string', list2, 'Max', length(list2), 'Min', 1)
popupmenuMap_Callback(handles.popupmenuMap, eventdata, handles)

% --- Executes on button press in pushbuttonAllState.
function pushbuttonAllState_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonAllState (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
list = get(handles.listboxStates, 'string');
set(handles.listboxPlotStates, 'string', list, 'Max', length(list), 'Min', 1)
popupmenuMap_Callback(handles.popupmenuMap, eventdata, handles)

% --- Executes on button press in pushbuttonClearState.
function pushbuttonClearState_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonClearState (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
list = {'None'};
set(handles.listboxPlotStates, 'string', list, 'Max', length(list), 'Min', 1)


% --- Executes on slider movement.
function sliderMarketPen_Callback(hObject, eventdata, handles)
% hObject    handle to sliderMarketPen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global result
result.Econ.MarketPen = round(get(handles.sliderMarketPen,'Value'));
set(handles.editMarketPen,'string',result.Econ.MarketPen)
SupplyDemand(handles)

% --- Executes during object creation, after setting all properties.
function sliderMarketPen_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderMarketPen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function editMarketPen_Callback(hObject, eventdata, handles)
% hObject    handle to editMarketPen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global result
result.Econ.MarketPen = str2double(get(handles.editMarketPen,'string'));
SupplyDemand(handles)

% --- Executes during object creation, after setting all properties.
function editMarketPen_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editMarketPen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function SupplyDemand(handles)
global result
buildNames = get(handles.listboxBuildings, 'string');
stateNames =  get(handles.listboxStates, 'string');
list1 = get(handles.listboxPlotBuild, 'string');
list2 = get(handles.listboxPlotStates, 'string');

load USbuildData
FCCost = zeros(50*16,1);
builds = zeros(50*16,1);
sizes = zeros(50*16,1);
for i = 1:1:length(list1)
    j = find(strcmp(list1(i),buildNames));
    numStates = length(list2);
    for k = 1:1:numStates
        state = find(strcmp(list2(k),stateNames));
        FCCost(k+numStates*(j-1),1) = result.FCevenCost(state,j);
        builds(k+numStates*(j-1),1) = buildData(state,j);
        sizes(k+numStates*(j-1),1) = result.FCsize(state,j);
    end
end

[DemandY, I] = sort(FCCost,1,'descend');

DemandX = zeros(length(FCCost),1);
DemandX(1) = builds(I(1))*sizes(I(1))/1000*result.Econ.MarketPen/100;
for ind = 2:1:length(FCCost)
    DemandX(ind) = DemandX(ind-1)+builds(I(ind))*sizes(I(ind))/1000*result.Econ.MarketPen/100;
end
axes(handles.axes2);
hold off
plot(DemandX,DemandY,'r') %convert to MW
hold on
plot(result.Econ.ProdCosts(:,1)/1000,result.Econ.ProdCosts(:,2))
axis([0 max(DemandX) 0 max(DemandY)])
xlabel('Annual Installed Capacity (MW)')
ylabel('Installed FC cost ($/kW)')
legend('Demand','Manufacture Cost','Location','NorthEast')


% --- Executes on button press in pushbuttonExport.
function pushbuttonExport_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonExport (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global result
ExportResult('NationalSurvey', result)