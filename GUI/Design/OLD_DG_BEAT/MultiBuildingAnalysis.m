function varargout = MultiBuildingAnalysis(varargin)
% MULTIBUILDINGANALYSIS M-file for MultiBuildingAnalysis.fig
%      MULTIBUILDINGANALYSIS, by itself, creates a new MULTIBUILDINGANALYSIS or raises the existing
%      singleton*.
%
%      H = MULTIBUILDINGANALYSIS returns the handle to a new MULTIBUILDINGANALYSIS or the handle to
%      the existing singleton*.
%
%      MULTIBUILDINGANALYSIS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MULTIBUILDINGANALYSIS.M with the given input arguments.
%
%      MULTIBUILDINGANALYSIS('Property','Value',...) creates a new MULTIBUILDINGANALYSIS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MultiBuildingAnalysis_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MultiBuildingAnalysis_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MultiBuildingAnalysis

% Last Modified by GUIDE v2.5 17-Mar-2014 17:21:22

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MultiBuildingAnalysis_OpeningFcn, ...
                   'gui_OutputFcn',  @MultiBuildingAnalysis_OutputFcn, ...
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


% --- Executes just before MultiBuildingAnalysis is made visible.
function MultiBuildingAnalysis_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MultiBuildingAnalysis (see VARARGIN)

% Choose default command line output for MultiBuildingAnalysis
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes MultiBuildingAnalysis wait for user response (see UIRESUME)
% uiwait(handles.figure1); 
global Project Model_dir
stateName = {'Alabama';'Alaska';'Arizona';'Arkansas';'California';'Colorado';'Connecticut';'Delaware';'Florida';'Georgia';
             'Hawaii';'Idaho';'Illinois';'Indiana';'Iowa';'Kansas';'Kentucky';'Louisiana';'Maine';'Maryland';
             'Massachusetts';'Michigan';'Minnesota';'Mississippi';'Missouri';'Montana';'Nebraska';'Nevada';'NewHampshire';'NewJersey';
             'NewMexico';'NewYork';'NorthCarolina';'NorthDakota';'Ohio';'Oklahoma';'Oregon';'Pennsylvania';'RhodeIsland';'SouthCarolina';
             'SouthDakota';'Tennessee';'Texas';'Utah';'Vermont';'Virginia';'Washington';'WestVirginia';'Wisconsin';'Wyoming';};
state = find(strcmp(Project.State,stateName));
set(handles.popupmenuState, 'string', stateName, 'value',state)
set(handles.editSellBack,'string',num2str(Project.Utilities.Grid.SellBackRate))
if ~isfield(Project.Utilities.Grid,'SellBackPerc')
    Project.Utilities.Grid.SellBackPerc =100;
end
set(handles.editSellbackPerc,'string',num2str(Project.Utilities.Grid.SellBackPerc))
if Project.Utilities.Grid.SellBackRate == 0
    set(handles.uipanelSellBack,'SelectedObject',handles.radiobuttonNoSellBack)
elseif Project.Utilities.Grid.SellBackRate == -1
    set(handles.uipanelSellBack,'SelectedObject',handles.radiobuttonReverseMeter)
else
    set(handles.uipanelSellBack,'SelectedObject',handles.radiobuttonFixedSellBack)
end
set(handles.uipanelSysSize,'SelectedObject',handles.radiobuttonOptSize)
set(handles.uipanelAnalysisType,'SelectedObject',handles.radiobuttonAllBuild)
set(handles.uipanelGridCosts,'SelectedObject',handles.radiobuttonScaleStateAvg)
if max(Project.Building.DemandC)>0
    set(handles.uipanelCHP,'SelectedObject',handles.radiobuttonDistHeatCool)
else set(handles.uipanelCHP,'SelectedObject',handles.radiobuttonNoDistHeatCool)
end



BuildType = {'Restaurant: full-service (sit down)';'Restaurant: quick-service (fast food)';'School: primary school';'School: secondary school';'Office: large office';'Office: medium office';'Office: small office';'Mid-rise apartment building';'Hospitality: large hotel';'Hospitality: small hotel/motel';'Health care: large hospital';'Health care: outpatient facility';'Retail: big-box, standalone retail store';'Retail: retail store located in a strip mall';'Retail: supermarket';'Unrefrigerated warehouse';};
Climate = {'Miami (ASHRAE 1A)';'Houston (ASHRAE 2A)';'Phoenix (ASHRAE 2B)';'Atlanta (ASHRAE 3A)';'Las Vegas (ASHRAE 3B-Inland)';'Los Angeles (ASHRAE 3B-Coast)';'San Francisco (ASHRAE 3C)';'Baltimore (ASHRAE 4A)';'Albuquerque (ASHRAE 4B)';'Seattle (ASHRAE 4C)';'Chicago (ASHRAE 5A)';'Boulder (ASHRAE 5B)';'Minneapolis (ASHRAE 6A)';'Helena, MT (ASHRAE 6B)';'Duluth, MN (ASHRAE 7)';'Fairbanks, AK (ASHRAE 8)';};
climateName = {'_1A_'; '_2A_'; '_2B_'; '_3A_'; '_3B_'; '_3B-Coast_'; '_3C_'; '_4A_'; '_4B_'; '_4C_'; '_5A_'; '_5B_'; '_6A_'; '_6B_'; '_7_'; '_8_';};
buildTypeName = {'SDRest'; 'FFRest'; 'Sch-pri'; 'Sch-sec'; 'LgOff'; 'MdOff'; 'SmOff'; 'MRapt'; 'LgHotel'; 'SmHotel'; 'Hospital'; 'OutP'; 'Retail'; 'StMall'; 'SMarket'; 'ware';};
Vintage = {'2010 construction (ASHRAE 90.1-2010)';'2007 construction (ASHRAE 90.1-2007)';'2004 construction 90.1-2004';'“Post-1980” construction (ASHRAE 90.1-1989)';'“Pre-1980” construction';};
vintageName = {'New2010';'New2007';'New2004';'Post1980';'Pre1980';};
name = char(Project.Building.NameList{1});
ind = strfind(name,'_');
build = find(strcmp(name(1:ind(1)-1),buildTypeName));
set(handles.popupmenuBuilding, 'string', BuildType, 'value',build)
clim = find(strcmp(name(ind(1):ind(2)),climateName));
set(handles.popupmenuClimate, 'string', Climate, 'value',clim) 
vin = find(strcmp(name(ind(2)+1:end),vintageName));
set(handles.popupmenuVintage, 'string', Vintage, 'value',vin) 
griddir=fullfile(Model_dir, 'System Library','Grid');
files=dir(fullfile(griddir,'*.mat'));
list=strrep({files.name},'.mat','');
utility = char(Project.Utilities.Grid.Name);
uti = find(strcmp(utility,list));
set(handles.popupmenuUtility, 'string', list, 'value',uti) 
controldir=fullfile(Model_dir, 'System Library','Control');
files=dir(fullfile(controldir,'*.mat'));
list=strrep({files.name},'.mat','');
control = char(Project.Control.Name);
con = find(strcmp(control,list));
set(handles.popupmenuControl, 'string', list, 'value',con) 


    
% --- Outputs from this function are returned to the command line.
function varargout = MultiBuildingAnalysis_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbuttonCancel.
function pushbuttonCancel_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonCancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clear global fixedSize testBuild
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
if get(handles.uipanelSysSize,'SelectedObject') == handles.radiobuttonFixSize
    SizeStrategy = 1;
elseif get(handles.uipanelSysSize,'SelectedObject') == handles.radiobuttonOptSize
    SizeStrategy = 2;
elseif get(handles.uipanelSysSize,'SelectedObject') == handles.radiobuttonBestNPV
    SizeStrategy = 3;
elseif get(handles.uipanelSysSize,'SelectedObject') == handles.radiobuttonLowestCO2
    SizeStrategy = 4;
end
if get(handles.uipanelCHP,'SelectedObject') == handles.radiobuttonDistHeatCool
    DistHeat = 1;
    DistCool = 1;
else DistHeat = 0;
    DistCool = 0;
end
defaultClimate = get(handles.popupmenuClimate, 'value');
climateName = {'_1A_'; '_2A_'; '_2B_'; '_3A_'; '_3B_'; '_3B-Coast_'; '_3C_'; '_4A_'; '_4B_'; '_4C_'; '_5A_'; '_5B_'; '_6A_'; '_6B_'; '_7_'; '_8_';};
if get(handles.uipanelAnalysisType,'SelectedObject') == handles.radiobuttonAllBuild || get(handles.uipanelAnalysisType,'SelectedObject') == handles.radiobuttonStates
    climate = climateName(defaultClimate);
else climate = climateName;
end

defaultBuilding = get(handles.popupmenuBuilding, 'value');
buildName = {'SDRest'; 'FFRest'; 'Sch-pri'; 'Sch-sec'; 'LgOff'; 'MdOff'; 'SmOff'; 'MRapt'; 'LgHotel'; 'SmHotel'; 'Hospital'; 'OutP'; 'Retail'; 'StMall'; 'SMarket'; 'ware';};
if get(handles.uipanelAnalysisType,'SelectedObject') == handles.radiobuttonAllClimate || get(handles.uipanelAnalysisType,'SelectedObject') == handles.radiobuttonStates
    buildType = buildName(defaultBuilding);
else buildType = buildName;
end

defaultState = get(handles.popupmenuState, 'value');
stateName = get(handles.popupmenuState, 'string');
if get(handles.uipanelAnalysisType,'SelectedObject') == handles.radiobuttonStates
    state = stateName;
else state = stateName(defaultState);
end

defaultVintage = get(handles.popupmenuVintage, 'value');
vintageName = {'New2010';'New2007';'New2004';'Post1980';'Pre1980';};
vintage = char(vintageName(defaultVintage));

% load utility if default has changed
UtilityList = get(handles.popupmenuUtility, 'string'); 
defaultUtility = char(UtilityList(get(handles.popupmenuUtility, 'value')));
if ~strcmp(Project.Utilities.Grid.Name,defaultUtility)
    griddir=fullfile(Model_dir, 'System Library','Grid');
    load(fullfile(griddir,defaultUtility));
    Project.Utilities.Grid = component;
end

% load control if default has changed
ControlList = get(handles.popupmenuControl, 'string'); 
defaultControl = char(ControlList(get(handles.popupmenuControl, 'value')));
if ~strcmp(Project.Control.Name,defaultControl)
    controldir=fullfile(Model_dir, 'System Library','Control');
    load(fullfile(controldir,defaultControl));
    Project.Control = component;
end

user = 'Commercial';
A = get(handles.uipanelGridCosts,'SelectedObject');
switch get(A,'Tag')
    case 'radiobuttonScaleStateAvg'
        scale = 1;
    case 'radiobuttonUtilitySpecific'
        scale = 0;
end

A = get(handles.uipanelGridEmissions,'SelectedObject');
switch get(A,'Tag')
    case 'radiobuttonStateAvg'
        GridMix = 1;
    case 'radiobuttonCombustionOnly'
        GridMix = 2;
end

name = char(Project.Building.NameList{1});
ind = strfind(name,'_');

count = 0;
h=waitbar(0,'Running Analysis');
for j = 1:length(climate)
    for k = 1:length(buildType)
        for s = 1:1:length(state)
            if length(state)>1
                climate(j) = char(ClimateZone(state(s),0)); %change climate each time
            end
            if ~strcmp(climate(j),name(ind(1):ind(2))) || ~strcmp(buildType(k),name(1:ind(1)-1))
                load(fullfile(Model_dir,'System Library','Buildings',strcat(buildType(k), climate(j), vintage)));
                %%% edit exterior lighting profile
                AbsMin =min(component.DemandE);
                LightLoad = max(component.ExteriorLight);
                AddLight = (component.ExteriorLight==0).*(component.DemandE<=(AbsMin+LightLoad));
                component.DemandE = component.DemandE+AddLight*LightLoad;
                %%% 
                Project.Building.DemandE = component.DemandE;
                Project.Building.DemandC = DistCool*component.DemandC;
                Project.Building.DemandH = DistHeat*component.DemandH;
                Project.Building.NonDistH = (1-DistHeat)*component.DemandH;
                Project.Building.CoolingElectricalLoad =  DistCool*component.CoolingElectricalLoad;
                Project.Building.NameList = cellstr(char(strcat(buildType(k),climate(j),vintage))); 
            end
            count = count+1;
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
            RESULT = runAnalyses(Project,user,scale);
            
            multiBuild.BuildType(count) = buildType(k);
            multiBuild.ClimateType(count) = climate(j);
            multiBuild.State(count) = state(s);
            multiBuild.Dispatch(count) = RESULT.Dispatch;
            multiBuild.Baseline(count) = RESULT.Baseline;
            multiBuild.costOut(count) = RESULT.costOut;
            multiBuild.eOut(count) = RESULT.eOut;
            
            % Demand Characteristics (describe shape of demand profile
            multiBuild.LoadFactor(count) = mean(Project.Building.DemandE)/max(Project.Building.DemandE);
            multiBuild.LoadScatter(count) = std(Project.Building.DemandE)/mean(Project.Building.DemandE);
            multiBuild.LoadEtoC(count) = mean(Project.Building.DemandE)/mean(Project.Building.DemandC);
            multiBuild.LoadEtoH(count) = mean(Project.Building.DemandE)/mean(Project.Building.DemandH);
            multiBuild.LoadCtoH(count) = mean(Project.Building.DemandC)/mean(Project.Building.DemandH);
            
            [BaselineEmission, DispatchEmission] = EmissionsCalculated(RESULT,state(s),GridMix,1);
            multiBuild.Summary.NPCbaseline(count,1:5) = [RESULT.costOut.NPVbaselineDemCharges RESULT.costOut.NPVbaselineUseCharges  RESULT.costOut.NPVbaselineFuelCost    RESULT.costOut.NPVbaselineOandM RESULT.costOut.NPVbaselineFinance;];
            multiBuild.Summary.NPCdispatch(count,1:5) = [RESULT.costOut.NPVnewDemCharges RESULT.costOut.NPVnewUseCharges RESULT.costOut.NPVnewFuelCost  RESULT.costOut.NPVnewOandM RESULT.costOut.NPVnewFinance;];
            multiBuild.Summary.CO2baseline(count,1:3) = BaselineEmission(1,1:3);
            multiBuild.Summary.NOxbaseline(count,1:3) = BaselineEmission(2,1:3);
            multiBuild.Summary.SO2baseline(count,1:3) = BaselineEmission(3,1:3);
            multiBuild.Summary.CO2dispatch(count,1:3) = DispatchEmission(1,1:3);
            multiBuild.Summary.NOxdispatch(count,1:3) = DispatchEmission(2,1:3);
            multiBuild.Summary.SO2dispatch(count,1:3) = DispatchEmission(3,1:3);
            
            waitbar(count/(length(climate)*length(buildType)*length(state)),h,strcat('Running Analysis  ',num2str(count),'  of  ',num2str(length(climate)*length(buildType)*length(state))));
        end
    end
end
close(h)

FCModel_dir=fullfile(Model_dir, 'results');
name = strcat('multiBuild',get(handles.editFileName,'String'));
filename=fullfile(FCModel_dir,name);
save(filename,'multiBuild')
close(gcf)
%% Visualization
MultiBuildResult()


% --- Executes on selection change in popupmenuBuilding.
function popupmenuBuilding_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuBuilding (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function popupmenuBuilding_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuBuilding (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenuClimate.
function popupmenuClimate_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuClimate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object creation, after setting all properties.
function popupmenuClimate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuClimate (see GCBO)
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

function editSellBack_Callback(hObject, eventdata, handles)
% hObject    handle to editSellBack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Project
Project.Utilities.Grid.SellBackRate = str2double(get(hObject,'String'));
if Project.Utilities.Grid.SellBackRate == 0
    set(handles.uipanelSellBack,'SelectedObject',handles.radiobuttonNoSellBack)
elseif Project.Utilities.Grid.SellBackRate == -1
    set(handles.uipanelSellBack,'SelectedObject',handles.radiobuttonReverseMeter)
else
    set(handles.uipanelSellBack,'SelectedObject',handles.radiobuttonFixedSellBack)
end

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


% --- Executes on selection change in popupmenuVintage.
function popupmenuVintage_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuVintage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object creation, after setting all properties.
function popupmenuVintage_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuVintage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenuUtility.
function popupmenuUtility_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuUtility (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object creation, after setting all properties.
function popupmenuUtility_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuUtility (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenuControl.
function popupmenuControl_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuControl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object creation, after setting all properties.
function popupmenuControl_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuControl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

