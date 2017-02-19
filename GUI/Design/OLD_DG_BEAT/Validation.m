function Validation(varargin)
% VALIDATION MATLAB code for Validation.fig
%      VALIDATION, by itself, creates a new VALIDATION or raises the existing
%      singleton*.
%
%      H = VALIDATION returns the handle to a new VALIDATION or the handle to
%      the existing singleton*.
%
%      VALIDATION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VALIDATION.M with the given input arguments.
%
%      VALIDATION('Property','Value',...) creates a new VALIDATION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Validation_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Validation_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Validation

% Last Modified by GUIDE v2.5 08-Apr-2014 17:39:22

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Validation_OpeningFcn, ...
                   'gui_OutputFcn',  @Validation_OutputFcn, ...
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


% --- Executes just before Validation is made visible.
function Validation_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Validation (see VARARGIN)

% Choose default command line output for Validation
global Project COMPONENT Model_dir
COMPONENT = Project.Building;
handles.output = hObject;

plotlist={'Electric Demand'
        'Electric without HVAC'
        'Cooling Demand'
        'Heat Demand'};
    
set(handles.checkboxDistHeat,'value',1);    
set(handles.checkboxDistCool,'value',1); 
set(handles.popupmenuAxes,'string',plotlist,'value',1)

resultlist= {'NPC'; 'CO2'; 'NOx'; 'SO2';};  
set(handles.popupmenuResult,'string',resultlist,'value',1)

builddir=fullfile(Model_dir, 'System Library','Buildings','RealBuildingData'); %strrep(which('NREL_FCModel.m'),'\main\NREL_FCModel.m','\component library\Building\RealBuildingData');
files=dir(fullfile(builddir,'*.mat'));
list=strrep({files.name},'.mat','');
set(handles.popupmenuDataFile,'string',list,'value',1)
popupmenuDataFile_Callback(handles.popupmenuDataFile,eventdata,handles,0);


controldir=fullfile(Model_dir, 'System Library','Control');% strrep(which('NREL_FCModel.m'),'\main\NREL_FCModel.m','\component library\Control');
files=dir(fullfile(controldir,'*.mat'));
list=strrep({files.name},'.mat','');
set(handles.popupmenuControl,'string',list,'value',1)
popupmenuControl_Callback(handles.popupmenuControl,eventdata,handles);

set(handles.sliderZoom,'Min',1,'Max',4,'Value',1)
sliderZoom_Callback(handles.sliderZoom, eventdata, handles)
listboxExistingBuild_Callback(handles.listboxExistingBuild,eventdata,handles)

% Update handles structure
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = Validation_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in listboxExistingBuild.
function listboxExistingBuild_Callback(hObject, eventdata, handles)
% hObject    handle to listboxExistingBuild (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%build list for Setup Listbox
global Project COMPONENT

vintageName = {'New2010';'New2007';'New2004';'Post1980';'Pre1980';};
climateName = {'_1A_'; '_2A_'; '_2B_'; '_3A_'; '_3B_'; '_3B-Coast_'; '_3C_'; '_4A_'; '_4B_'; '_4C_'; '_5A_'; '_5B_'; '_6A_'; '_6B_'; '_7_'; '_8_';};
buildTypeName = {'SDRest'; 'FFRest'; 'Sch-pri'; 'Sch-sec'; 'LgOff'; 'MdOff'; 'SmOff'; 'MRapt'; 'LgHotel'; 'SmHotel'; 'Hospital'; 'OutP'; 'Retail'; 'StMall'; 'SMarket'; 'ware';};

BuildType = {'Restaurant: full-service (sit down)';'Restaurant: quick-service (fast food)';'School: primary school';'School: secondary school';'Office: large office';'Office: medium office';'Office: small office';'Mid-rise apartment building';'Hospitality: large hotel';'Hospitality: small hotel/motel';'Health care: large hospital';'Health care: outpatient facility';'Retail: big-box, standalone retail store';'Retail: retail store located in a strip mall';'Retail: supermarket';'Unrefrigerated warehouse';};
Climate = {'Miami (ASHRAE 1A)';'Houston (ASHRAE 2A)';'Phoenix (ASHRAE 2B)';'Atlanta (ASHRAE 3A)';'Las Vegas (ASHRAE 3B-Inland)';'Los Angeles (ASHRAE 3B-Coast)';'San Francisco (ASHRAE 3C)';'Baltimore (ASHRAE 4A)';'Albuquerque (ASHRAE 4B)';'Seattle (ASHRAE 4C)';'Chicago (ASHRAE 5A)';'Boulder (ASHRAE 5B)';'Minneapolis (ASHRAE 6A)';'Helena, MT (ASHRAE 6B)';'Duluth, MN (ASHRAE 7)';'Fairbanks, AK (ASHRAE 8)';};
Vintage = {'2010 construction (ASHRAE 90.1-2010)';'2007 construction (ASHRAE 90.1-2007)';'2004 construction 90.1-2004';'“Post-1980” construction (ASHRAE 90.1-1989)';'“Pre-1980” construction';};
set(handles.lbBuildingType, 'string', BuildType(:,1), 'Max', length(BuildType(:,1)), 'Min', 1)
set(handles.lbBuildingClimate, 'string', Climate(:,1), 'Max', length(Climate(:,1)), 'Min', 1)
set(handles.lbBuildingVintage, 'string', Vintage(:,1), 'Max', length(Vintage(:,1)), 'Min', 1)

    %build list of existing names
str={};
for i = 1:1:length(COMPONENT.NameList)
    name = char(COMPONENT.NameList{i});
    ind =strfind(name,'_');
    build = find(strcmp(name(1:ind(1)-1),buildTypeName));
    clim = find(strcmp(name(ind(1):ind(2)),climateName));
    vin = find(strcmp(name(ind(2)+1:end),vintageName));
    str(end+1) = {char(strcat(BuildType(build), Climate(clim), Vintage(vin)))};
end


if length(COMPONENT.NameList) == 1
    name = char(COMPONENT.NameList{1});
else name = 'Multiple';
end
COMPONENT.Name = name;
if ~isempty(COMPONENT.NameList)
    Project.Building = COMPONENT;
    Project.Result = runAnalyses(Project);
    popupmenuAxes_Callback(handles.popupmenuAxes,eventdata,handles)
    popupmenuResult_Callback(handles.popupmenuResult, eventdata, handles)    
end
set(handles.listboxExistingBuild,'string',char(str))
CHPsize = 0;
for i = 1:1:length(Project.System.CHP)
    CHPsize = CHPsize + Project.System.CHP(i).SysSize(1);
end
set(handles.editCHPsize,'string',CHPsize)


% --- Executes during object creation, after setting all properties.
function listboxExistingBuild_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listboxExistingBuild (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in lbBuildingType.
function lbBuildingType_Callback(hObject, eventdata, handles)
% hObject    handle to lbBuildingType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
popupmenuAxes_Callback(handles.popupmenuAxes,eventdata,handles)

% --- Executes during object creation, after setting all properties.
function lbBuildingType_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lbBuildingType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in lbBuildingClimate.
function lbBuildingClimate_Callback(hObject, eventdata, handles)
% hObject    handle to lbBuildingClimate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
popupmenuAxes_Callback(handles.popupmenuAxes,eventdata,handles)

% --- Executes during object creation, after setting all properties.
function lbBuildingClimate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lbBuildingClimate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in lbBuildingVintage.
function lbBuildingVintage_Callback(hObject, eventdata, handles)
% hObject    handle to lbBuildingVintage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
popupmenuAxes_Callback(handles.popupmenuAxes,eventdata,handles)

% --- Executes during object creation, after setting all properties.
function lbBuildingVintage_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lbBuildingVintage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pushbuttonDone.
function pushbuttonDone_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonDone (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Project COMPONENT
Project.Building = COMPONENT;
close(gcf)
uiresume


function popupmenuAxes_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuAxes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
val=get(hObject,'value');
str=get(hObject,'string');
currentstr=str{val};

global Project RealBuild Model_dir
Ts = 8760/length(Project.Building.DemandE);
axes(findobj(gcf,'tag','axes1'))
cla reset
set(gca,'tag','axes1')
hold all

result = Project.Result.eOut;

%load new building profile (this is the building selected but not yet added)
type=get(handles.lbBuildingType, 'value');
climate=get(handles.lbBuildingClimate, 'value');
vintage=get(handles.lbBuildingVintage, 'value');
type = type(1);
climate = climate(1);
vintage = vintage(1);
vintageName = {'New2010';'New2007';'New2004';'Post1980';'Pre1980';};
climateName = {'_1A_'; '_2A_'; '_2B_'; '_3A_'; '_3B_'; '_3B-Coast_'; '_3C_'; '_4A_'; '_4B_'; '_4C_'; '_5A_'; '_5B_'; '_6A_'; '_6B_'; '_7_'; '_8_';};
buildTypeName = {'SDRest'; 'FFRest'; 'Sch-pri'; 'Sch-sec'; 'LgOff'; 'MdOff'; 'SmOff'; 'MRapt'; 'LgHotel'; 'SmHotel'; 'Hospital'; 'OutP'; 'Retail'; 'StMall'; 'SMarket'; 'ware';};
name = char(strcat(buildTypeName(type),climateName(climate),vintageName(vintage)));
buildDir = fullfile(Model_dir, 'System Library' ,'Buildings'); % strrep(which('NREL_FCModel.m'),'\main\NREL_FCModel.m','\component library\Building');
load(strcat(buildDir,filesep,name))
newBuild = component;
%%%

Zoom = max(1,round(get(handles.sliderZoom,'Value')));
StartDate = max(1,floor(get(handles.sliderDate,'Value')));
month = [0 31 28 31 30 31 30 31 31 30 31 30 31];
monthDays = [0 31 59 90 120 151 181 212 243 273 304 334 365];
monthLabel = ['January  '; 'February '; 'March    '; 'April    '; 'May      '; 'June     '; 'July     '; 'August   '; 'September'; 'October  ';'November ' ;'December ';];
if Zoom == 1;
    day1 = 1;
    lastDay = 365;
    plotAxis = datenum([2014*ones(12,1) (1:12)' 1*ones(12,1) 0*ones(12,1) 0*ones(12,1) 0*ones(12,1)]);
    xlabel('Date')
elseif Zoom == 2;
    day1 = 1+sum(month(1:StartDate));
    lastDay = sum(month(1:StartDate+1));
    days = month(StartDate+1);
    plotAxis = datenum([2014*ones(days,1) StartDate*ones(days,1) (1:days)' 0*ones(days,1) 0*ones(days,1) 0*ones(days,1)]);
	xlabel(monthLabel(StartDate,:))
elseif Zoom == 3;
    day1 = 1+7*(StartDate-1);
    lastDay = 1+7*StartDate;
    StartMonth = find(monthDays>day1,1) - 1;
    StartDay = day1-monthDays(StartMonth);
    plotAxis = datenum([2014*ones(8,1) StartMonth*ones(8,1)  (StartDay:StartDay+7)'  0*ones(8,1)  0*ones(8,1) 0*ones(8,1)]);
    if StartDay+6<=month(StartMonth+1)
        xlabel(monthLabel(StartMonth,:))
    else 
        xlabel(strcat(monthLabel(StartMonth,:), ' / ', monthLabel(StartMonth+1,:)))
    end
elseif Zoom == 4;
    day1 = StartDate;
    lastDay = StartDate;
    StartMonth = find(monthDays>day1,1) - 1;
    StartDay = day1-monthDays(StartMonth);
    plotAxis = datenum([2014*ones(24,1)  StartMonth*ones(24,1)  StartDay*ones(24,1) (1:24)' 0*ones(24,1) 0*ones(24,1)]);
    xlabel(strcat(['Hours of ', monthLabel(StartMonth,:), num2str(StartDay)]))
end
X = datenum(2014,1,1,0,0,0)+((day1-1)+(Ts/24):(Ts/24):lastDay);
switch currentstr
    case 'Heat Demand'
        Y1 = RealBuild.DemandH(1+24/Ts*(day1-1):24/Ts*lastDay);
        Y2 = Project.Building.DemandH(1+24/Ts*(day1-1):24/Ts*lastDay);
        Y3 = newBuild.DemandH(1+24/Ts*(day1-1):24/Ts*lastDay);
    case 'Electric Demand'
        Y1 = RealBuild.DemandE(1+24/Ts*(day1-1):24/Ts*lastDay);
        Y2 = Project.Building.DemandE(1+24/Ts*(day1-1):24/Ts*lastDay);
        Y3 = newBuild.DemandE(1+24/Ts*(day1-1):24/Ts*lastDay);
    case 'Cooling Demand'
        Y1 = RealBuild.DemandC(1+24/Ts*(day1-1):24/Ts*lastDay);
        Y2 = Project.Building.DemandC(1+24/Ts*(day1-1):24/Ts*lastDay);
        Y3 = newBuild.DemandC(1+24/Ts*(day1-1):24/Ts*lastDay);
    case 'Electric without HVAC'
        Y1 = RealBuild.DemandE(1+24/Ts*(day1-1):24/Ts*lastDay) - RealBuild.CoolingElectricalLoad(1+24/Ts*(day1-1):24/Ts*lastDay);
        Y2 = Project.Building.DemandE(1+24/Ts*(day1-1):24/Ts*lastDay)-Project.Building.CoolingElectricalLoad(1+24/Ts*(day1-1):24/Ts*lastDay);
        Y3 = newBuild.DemandE(1+24/Ts*(day1-1):24/Ts*lastDay)-newBuild.CoolingElectricalLoad(1+24/Ts*(day1-1):24/Ts*lastDay);      
end
if get(handles.checkboxReal,'Value')==1
    plot(X,Y1);
end
if get(handles.checkboxSim,'Value')==1
    plot(X,Y2);
end
if get(handles.checkboxSelected,'Value')==1
    plot(X,Y3);
end
ylabel('Demand') 
set(gca,'xtick',plotAxis) 
if Zoom ==1
    datetick('x','mmmdd','keepticks')
elseif Zoom ==2
    datetick('x','dd','keepticks')
elseif Zoom ==3
    datetick('x','dd','keepticks')
elseif Zoom ==4
    datetick('x','HH','keepticks')
end
str = get(hObject,'string');
elec = strcmp(str(get(handles.popupmenuAxes,'value')),'Electric Demand') || strcmp(str(get(handles.popupmenuAxes,'value')),'Electric without HVAC');
if elec
    set(handles.checkboxGeneration,'enable','on')
else set(handles.checkboxGeneration,'enable','off')
end
if isfield(Project.System,'TES')
    set(handles.checkboxTESshift,'enable','on')
else set(handles.checkboxTESshift,'enable','off')
end
if strcmp(str(get(hObject,'value')),'Electric Demand') || strcmp(str(get(hObject,'value')),'Electric without HVAC')
    if get(handles.checkboxGeneration,'Value')==1
        genPow = result.GenPower(1+24/Ts*(day1-1):24/Ts*lastDay);
        hold on
        plot(X,genPow,'r');
        hold off
    end
    if get(handles.checkboxTESshift,'Value')==1 && isfield(Project.System,'TES')
        hold on
        plot(X,result.DemandEshift(1+24/Ts*(day1-1):24/Ts*lastDay),'g');
        hold off
    end
end

% --- Executes on slider movement.
function sliderZoom_Callback(hObject, eventdata, handles)
% hObject    handle to sliderZoom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Zoom = round(get(hObject,'Value'));
if Zoom == 1
    set(handles.sliderDate,'Min',1)
    set(handles.sliderDate,'Max',1.1)
elseif Zoom == 2
    set(handles.sliderDate,'Min',1)
    set(handles.sliderDate,'Max',12)
elseif Zoom == 3
    set(handles.sliderDate,'Min',1)
    set(handles.sliderDate,'Max',52)
elseif Zoom == 4
    set(handles.sliderDate,'Min',1)
    set(handles.sliderDate,'Max',365)
end
set(handles.sliderDate,'Value',1)
popupmenuAxes_Callback(handles.popupmenuAxes,eventdata,handles)

% --- Executes during object creation, after setting all properties.
function popupmenuAxes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuAxes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function sliderZoom_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderZoom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function sliderDate_Callback(hObject, eventdata, handles)
% hObject    handle to sliderDate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
popupmenuAxes_Callback(handles.popupmenuAxes,eventdata,handles)

% --- Executes during object creation, after setting all properties.
function sliderDate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderDate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes on button press in checkboxGeneration.
function checkboxGeneration_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxGeneration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
popupmenuAxes_Callback(handles.popupmenuAxes,eventdata,handles)


% --- Executes on button press in checkboxTESshift.
function checkboxTESshift_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxTESshift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
popupmenuAxes_Callback(handles.popupmenuAxes,eventdata,handles)


% --- Executes on button press in pushbuttonAdd.
function pushbuttonAdd_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonAdd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global COMPONENT Model_dir 
type=get(handles.lbBuildingType, 'value');
climate=get(handles.lbBuildingClimate, 'value');
vintage=get(handles.lbBuildingVintage, 'value');

vintageName = {'New2010';'New2007';'New2004';'Post1980';'Pre1980';};
climateName = {'_1A_'; '_2A_'; '_2B_'; '_3A_'; '_3B_'; '_3B-Coast_'; '_3C_'; '_4A_'; '_4B_'; '_4C_'; '_5A_'; '_5B_'; '_6A_'; '_6B_'; '_7_'; '_8_';};
buildTypeName = {'SDRest'; 'FFRest'; 'Sch-pri'; 'Sch-sec'; 'LgOff'; 'MdOff'; 'SmOff'; 'MRapt'; 'LgHotel'; 'SmHotel'; 'Hospital'; 'OutP'; 'Retail'; 'StMall'; 'SMarket'; 'ware';};
buildDir =  fullfile(Model_dir,'System Library','Buildings'); %strrep(which('NREL_FCModel.m'),'\main\NREL_FCModel.m','\component library\Building');
for i = 1:1:length(type)
    for j = 1:1:length(climate)
        for k = 1:1:length(vintage)
            Ti = type(i);
            Cj = climate(j);
            Vk = vintage(k);
            name = char(strcat(buildTypeName(Ti),climateName(Cj),vintageName(Vk)));
            load(strcat(buildDir,filesep,name))
            COMPONENT.DistHeat(end+1) = get(handles.checkboxDistHeat,'value');
            COMPONENT.DistCool(end+1) = get(handles.checkboxDistCool,'value');
            COMPONENT.DemandE = COMPONENT.DemandE + component.DemandE;
            COMPONENT.DemandC = COMPONENT.DemandC + COMPONENT.DistCool(end)*component.DemandC;
            COMPONENT.DemandH = COMPONENT.DemandH + COMPONENT.DistHeat(end)*component.DemandH;
            COMPONENT.NonDistH = COMPONENT.NonDistH + (1-COMPONENT.DistHeat(end))*component.DemandH;
            COMPONENT.CoolingElectricalLoad = COMPONENT.CoolingElectricalLoad + COMPONENT.DistCool(end)*component.CoolingElectricalLoad;
            COMPONENT.NameList(end+1) = cellstr(name);  
        end
    end
end
listboxExistingBuild_Callback(handles.listboxExistingBuild, eventdata, handles)

% --- Executes on button press in pushbuttonRemove.
function pushbuttonRemove_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonRemove (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global COMPONENT Model_dir
val = get(handles.listboxExistingBuild,'value');
buildDir =  fullfile(Model_dir,'System Library','Buildings'); %strrep(which('NREL_FCModel.m'),'\main\NREL_FCModel.m','\component library\Building');
for j = 1:1:length(val)
    name = char(COMPONENT.NameList(val(j)));
    load(strcat(buildDir,filesep,name))
    COMPONENT.DemandE = COMPONENT.DemandE - component.DemandE;
    COMPONENT.DemandC = COMPONENT.DemandC - COMPONENT.DistCool(val(j))*component.DemandC;
    COMPONENT.DemandH = COMPONENT.DemandH - COMPONENT.DistHeat(val(j))*component.DemandH;
    COMPONENT.NonDistH = COMPONENT.NonDistH - (1-COMPONENT.DistHeat(val(j)))*component.DemandH;
    COMPONENT.CoolingElectricalLoad = COMPONENT.CoolingElectricalLoad - COMPONENT.DistCool(val(j))*component.CoolingElectricalLoad;
end
count = 0;
AllNames = COMPONENT.NameList;
DistHeat = COMPONENT.DistHeat;
DistCool = COMPONENT.DistCool;
COMPONENT.NameList = {};
COMPONENT.DistHeat = [];
COMPONENT.DistCool = [];
for i = 1:1:length(AllNames)
    if i ~=val
        count = count+1;
        COMPONENT.NameList(count) = AllNames(i);
        COMPONENT.DistHeat(count) = DistHeat(i);
        COMPONENT.DistCool(count) = DistCool(i);
    end
end
set(handles.listboxExistingBuild,'Value',1)
listboxExistingBuild_Callback(handles.listboxExistingBuild, eventdata, handles)

% --- Executes on button press in pushbuttonClear.
function pushbuttonClear_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonClear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global COMPONENT
Ts = 8760/length(COMPONENT.DemandE);
COMPONENT.NameList = {};
COMPONENT.DistHeat = [];
COMPONENT.DistCool = [];
COMPONENT.DemandE = zeros(8760/Ts,1);
COMPONENT.DemandC = zeros(8760/Ts,1);
COMPONENT.DemandH = zeros(8760/Ts,1);
COMPONENT.NonDistH = zeros(8760/Ts,1);
COMPONENT.CoolingElectricalLoad = zeros(8760/Ts,1);
set(handles.listboxExistingBuild,'Value',1)
listboxExistingBuild_Callback(handles.listboxExistingBuild, eventdata, handles)

% --- Executes on button press in pushbuttonSwitch.
function pushbuttonSwitch_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonSwitch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
pushbuttonRemove_Callback(hObject, eventdata, handles)
pushbuttonAdd_Callback(hObject, eventdata, handles)

% --- Executes on button press in pushbuttonScale.
function pushbuttonScale_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonScale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global COMPONENT RealBuild
avgEreal = mean(RealBuild.DemandE);
avgEsim = mean(COMPONENT.DemandE);
ScaleRatio = avgEreal/avgEsim;
COMPONENT.DemandE = COMPONENT.DemandE*ScaleRatio;
COMPONENT.DemandC = COMPONENT.DemandC*ScaleRatio;
COMPONENT.DemandH = COMPONENT.DemandH*ScaleRatio;
COMPONENT.CoolingElectricalLoad = COMPONENT.CoolingElectricalLoad*ScaleRatio;
listboxExistingBuild_Callback(handles.listboxExistingBuild, eventdata, handles)

% --- Executes on button press in checkboxDistHeat.
function checkboxDistHeat_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxDistHeat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on button press in checkboxDistCool.
function checkboxDistCool_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxDistCool (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on selection change in popupmenuDataFile.
function popupmenuDataFile_Callback(hObject, eventdata, handles,skipLoad)
% hObject    handle to popupmenuDataFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global RealBuild Model_dir
str=get(hObject,'string');
val=get(hObject,'value');
if nargin==3 || skipLoad ==0
    dirToLoad=fullfile(Model_dir,'System Library','Building', 'RealBuildingData');
    load(fullfile(dirToLoad,str{val}))
end
RealBuild = A;
RealBuild.NonDistH = zeros(length(A.DemandE(:,1)),1);
RealBuild.CHPtemp = 80;
popupmenuAxes_Callback(handles.popupmenuAxes,eventdata,handles)
popupmenuResult_Callback(handles.popupmenuResult, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function popupmenuDataFile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuDataFile (see GCBO)
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
global Project Model_dir
str=get(hObject,'string');
val=get(hObject,'value');
dirToLoad=fullfile(Model_dir,'System Library','Control');
load(fullfile(dirToLoad,str{val}));
Project.Control = component;
listboxExistingBuild_Callback(handles.listboxExistingBuild, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function popupmenuControl_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuControl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenuResult.
function popupmenuResult_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuResult (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Project RealBuild
val=get(hObject,'value');
str=get(hObject,'string');
currentstr=str{val};
axes(handles.axes2);
if isfield(Project,'Renewable')
    renew = Project.Renewable;
else renew = 1;
end

%Baseline = Real Building Data
DispatchRealBuild=dispatchSystem(Project.System, RealBuild,Project.Utilities.Grid,renew,Project.Control,'California');

%Baseline
RESULT.Baseline.Elec = zeros(8760,1);
steps = length(RealBuild.DemandE);
Ts = 8760/steps;
for i = 1:1:8760
    RESULT.Baseline.Elec(i) = sum(RealBuild.DemandE(1+1/Ts*(i-1):1/Ts*i)*Ts);
    RESULT.Baseline.Heat(i) = sum(RealBuild.DemandH(1+1/Ts*(i-1):1/Ts*i)*Ts);
end
monthDays = [0 31 59 90 120 151 181 212 243 273 304 334 365];
for i = 1:1:length(monthDays)-1
    RESULT.Baseline.Fuel(i) = sum(RealBuild.DemandH(24/Ts*monthDays(i)+1:24/Ts*monthDays(i+1)))*Ts;
    RESULT.Dispatch.Fuel(i) = sum(DispatchRealBuild.eOut.Fuel(24/Ts*monthDays(i)+1:24/Ts*monthDays(i+1)));
end
%Dispatch
RESULT.Dispatch.Elec = DispatchRealBuild.eOut.Grid_purchases_hourly_kWh;
RESULT.Dispatch.ElecTotProd = DispatchRealBuild.eOut.Total_Electricity_produced_kWh;

RESULT.eOut = DispatchRealBuild.eOut;
RESULT.Dispatch.SysSize = DispatchRealBuild.eOut.SysSize;
RESULT.Dispatch.ChillerSize(1) = DispatchRealBuild.eOut.ChillerSize(1);
RESULT.Dispatch.ChillerSize(2) = DispatchRealBuild.eOut.ChillerSize(2);
RESULT.Dispatch.TESsize = DispatchRealBuild.eOut.TESsize;
RESULT.Dispatch.BatterySize = DispatchRealBuild.eOut.BatterySize;
RESULT.Dispatch.RenewableSize(1) = DispatchRealBuild.eOut.RenewableSize(1);
RESULT.Dispatch.RenewableSize(2) = DispatchRealBuild.eOut.RenewableSize(2); 
costRealBuild=FinancialCalcs2(RESULT.Baseline,RESULT.Dispatch,Project.Utilities.Grid,Project.Economic);

switch currentstr
    case 'NPC'
        bar([costRealBuild.NPVbaselineDemCharges costRealBuild.NPVbaselineUseCharges costRealBuild.NPVbaselineFuelCost  costRealBuild.NPVbaselineOandM costRealBuild.NPVbaselineFinance;...
             costRealBuild.NPVnewDemCharges costRealBuild.NPVnewUseCharges costRealBuild.NPVnewFuelCost  costRealBuild.NPVnewOandM costRealBuild.NPVnewFinance],'stacked')
        ylabel('Fuel, Grid, and O&M/Finance Costs ($)')
        set(gca,'XTickLabel',{'Baseline','System'})
        xlim([0.5 2.5])
        legend('Demand', 'Grid' ,'Fuel', 'O&M','Financing','Location','NorthEastOutside')
	case 'CO2'
        [BaselineEmission DispatchEmission] = EmissionsCalculated(RESULT,'CA',2,1);
        bar([BaselineEmission(1,1:3);DispatchEmission(1,1:3);],'stacked')
        ylabel('Emissions of CO2 (tons)')
        set(gca,'XTickLabel',{'Baseline','System'})
        xlim([0.5 2.5])
        legend('Grid', 'CHP' ,'Boiler', 'Location','NorthEastOutside')
    case 'NOx'
        [BaselineEmission DispatchEmission] = EmissionsCalculated(RESULT,'CA',2,1);
        bar([BaselineEmission(2,1:3);DispatchEmission(2,1:3);],'stacked')
        ylabel('Emissions of  NOx (lbs)')
        set(gca,'XTickLabel',{'Baseline','System'})
        xlim([0.5 2.5])
        legend('Grid', 'CHP' ,'Boiler', 'Location','NorthEastOutside')
    case 'SO2'
        [BaselineEmission DispatchEmission] = EmissionsCalculated(RESULT,'CA',2,1);
        bar([BaselineEmission(3,1:3);DispatchEmission(3,1:3);],'stacked')
        ylabel('Emissions of SO2 (lbs)')
        set(gca,'XTickLabel',{'Baseline','System'})
        xlim([0.5 2.5])
        legend('Grid', 'CHP' ,'Boiler', 'Location','NorthEastOutside')
end



% --- Executes during object creation, after setting all properties.
function popupmenuResult_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuResult (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkboxReal.
function checkboxReal_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxReal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
popupmenuAxes_Callback(handles.popupmenuAxes,eventdata,handles)

% --- Executes on button press in checkboxSim.
function checkboxSim_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxSim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
popupmenuAxes_Callback(handles.popupmenuAxes,eventdata,handles)

% --- Executes on button press in checkboxSelected.
function checkboxSelected_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxSelected (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
popupmenuAxes_Callback(handles.popupmenuAxes,eventdata,handles)

function editCHPsize_Callback(hObject, eventdata, handles)
% hObject    handle to editCHPsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Project
for i = 1:1:length(Project.System.CHP)
    oldSize(i) = Project.System.CHP(i).SysSize(1);
end
AllSize = sum(oldSize);
NewSize = str2double(get(hObject,'String'));
for i = 1:1:length(Project.System.CHP)
    Project.System.CHP(i).SysSize(1) = oldSize(i)*(NewSize/AllSize);
end
listboxExistingBuild_Callback(handles.listboxExistingBuild, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function editCHPsize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editCHPsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pushbuttonOptSize.
function pushbuttonOptSize_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonOptSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ScaleSystemComponents(100,100,100,120,1,2,1,2,2)
listboxExistingBuild_Callback(handles.listboxExistingBuild, eventdata, handles)


% --- Executes on button press in pushbuttonCostSize.
function pushbuttonCostSize_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonCostSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Project
if isfield(Project.System,'TES')
    autoSizing('cost','TES')
end
autoSizing('cost','CHP')
if isfield(Project.System,'Battery')
    autoSizing('cost','Battery',3)
end
listboxExistingBuild_Callback(handles.listboxExistingBuild, eventdata, handles)


% --- Executes on button press in pushbuttonResult.
function pushbuttonResult_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonResult (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global RealBuild Project Model_dir
Project.Building = RealBuild;
Project.Building.NameList = RealBuild.Name;
Project.Building.DistHeat(1) = get(handles.checkboxDistHeat,'value');
Project.Building.DistCool(1) = get(handles.checkboxDistCool,'value');
projdir=fullfile(Model_dir, 'project');
Project.Result = runAnalyses(Project);
save(fullfile(projdir,strcat(char(RealBuild.Name),'.mat')),'Project')
files=dir(fullfile(projdir,'*.mat'));
list=strrep({files.name},'.mat','');
currentProject = find(strcmp(Project.Name,list));
close(gcf)
uiresume
ViewResults2(currentProject)
