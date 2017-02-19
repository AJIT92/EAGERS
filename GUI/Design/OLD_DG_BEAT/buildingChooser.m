function buildingChooser(varargin)
% BUILDINGCHOOSER MATLAB code for buildingChooser.fig
%      BUILDINGCHOOSER, by itself, creates a new BUILDINGCHOOSER or raises the existing
%      singleton*.
%
%      H = BUILDINGCHOOSER returns the handle to a new BUILDINGCHOOSER or the handle to
%      the existing singleton*.
%
%      BUILDINGCHOOSER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BUILDINGCHOOSER.M with the given input arguments.
%
%      BUILDINGCHOOSER('Property','Value',...) creates a new BUILDINGCHOOSER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before buildingChooser_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to buildingChooser_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help buildingChooser

% Last Modified by GUIDE v2.5 27-Jan-2014 16:38:30

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @buildingChooser_OpeningFcn, ...
                   'gui_OutputFcn',  @buildingChooser_OutputFcn, ...
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


% --- Executes just before buildingChooser is made visible.
function buildingChooser_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to buildingChooser (see VARARGIN)

% Choose default command line output for buildingChooser
global Project COMPONENT
handles.output = hObject;
COMPONENT = Project.Building;

plotlist={'Electric Demand'
        'Electric without HVAC'
        'Electric Demand by Day'
        'Cooling Demand'
        'Heat Demand'
        'New Elec Demand'
        'New Elec w/o HVAC'
        'New Cooling Demand'};
    
set(handles.checkboxDistHeat,'value',1);    
set(handles.checkboxDistCool,'value',1); 
set(handles.popupmenuAxes,'string',plotlist,'value',1)
set(handles.sliderZoom,'Min',1,'Max',4,'Value',1,'SliderStep',[1/3,1/3])
sliderZoom_Callback(handles.sliderZoom, eventdata, handles)

BuildType = {'Restaurant: full-service (sit down)';'Restaurant: quick-service (fast food)';'School: primary school';'School: secondary school';'Office: large office';'Office: medium office';'Office: small office';'Mid-rise apartment building';'Hospitality: large hotel';'Hospitality: small hotel/motel';'Health care: large hospital';'Health care: outpatient facility';'Retail: big-box, standalone retail store';'Retail: retail store located in a strip mall';'Retail: supermarket';'Unrefrigerated warehouse';};
Climate = {'Miami (ASHRAE 1A)';'Houston (ASHRAE 2A)';'Phoenix (ASHRAE 2B)';'Atlanta (ASHRAE 3A)';'Las Vegas (ASHRAE 3B-Inland)';'Los Angeles (ASHRAE 3B-Coast)';'San Francisco (ASHRAE 3C)';'Baltimore (ASHRAE 4A)';'Albuquerque (ASHRAE 4B)';'Seattle (ASHRAE 4C)';'Chicago (ASHRAE 5A)';'Boulder (ASHRAE 5B)';'Minneapolis (ASHRAE 6A)';'Helena, MT (ASHRAE 6B)';'Duluth, MN (ASHRAE 7)';'Fairbanks, AK (ASHRAE 8)';};
Vintage = {'2010 construction (ASHRAE 90.1-2010)';'2007 construction (ASHRAE 90.1-2007)';'2004 construction 90.1-2004';'“Post-1980” construction (ASHRAE 90.1-1989)';'“Pre-1980” construction';};
set(handles.lbBuildingType, 'string', BuildType(:,1), 'Max', length(BuildType(:,1)), 'Min', 1)
set(handles.lbBuildingClimate, 'string', Climate(:,1), 'Max', length(Climate(:,1)), 'Min', 1)
set(handles.lbBuildingVintage, 'string', Vintage(:,1), 'Max', length(Vintage(:,1)), 'Min', 1)

for i = 1:1:length(COMPONENT.DistHeat)
    if ~exist('COMPONENT.CHPtemp(i)')
        COMPONENT.CHPtemp(i) = 80;
        CHPtemp(i) = 80;
    else
        CHPtemp(i) = COMPONENT.CHPtemp(i);
    end
end
CHPtempMin = min(CHPtemp);
set(handles.editCHPtemp,'string',num2str(CHPtempMin));

listboxExistingBuild_Callback(handles.listboxExistingBuild,eventdata,handles)

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes buildingChooser wait for user response (see UIRESUME)


% --- Outputs from this function are returned to the command line.
function varargout = buildingChooser_OutputFcn(hObject, eventdata, handles) 
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
global COMPONENT

vintageName = {'New2010';'New2007';'New2004';'Post1980';'Pre1980';};
climateName = {'_1A_'; '_2A_'; '_2B_'; '_3A_'; '_3B_'; '_3B-Coast_'; '_3C_'; '_4A_'; '_4B_'; '_4C_'; '_5A_'; '_5B_'; '_6A_'; '_6B_'; '_7_'; '_8_';};
buildTypeName = {'SDRest'; 'FFRest'; 'Sch-pri'; 'Sch-sec'; 'LgOff'; 'MdOff'; 'SmOff'; 'MRapt'; 'LgHotel'; 'SmHotel'; 'Hospital'; 'OutP'; 'Retail'; 'StMall'; 'SMarket'; 'ware';};

BuildType = get(handles.lbBuildingType, 'string');
Climate = get(handles.lbBuildingClimate, 'string');
Vintage = get(handles.lbBuildingVintage, 'string');

    %build list of existing names
str={};
if ~isfield(COMPONENT,'NameList')
    str = char(COMPONENT.Name);
else
    for i = 1:1:length(COMPONENT.NameList)
        name = char(COMPONENT.NameList{i});
        ind =strfind(name,'_');
        build = find(strcmp(name(1:ind(1)-1),buildTypeName));
        clim = find(strcmp(name(ind(1):ind(2)),climateName));
        vin = find(strcmp(name(ind(2)+1:end),vintageName));
        str(end+1) = {char(strcat(BuildType(build), Climate(clim), Vintage(vin)))};
    end
end
if ~isfield(COMPONENT,'NameList')
    name = char(COMPONENT.Name);
elseif length(COMPONENT.NameList) == 1
    name = char(COMPONENT.NameList{1});
else name = 'Multiple';
end
COMPONENT.Name = name;
if isfield(COMPONENT,'NameList') && ~isempty(COMPONENT.NameList)
    popupmenuAxes_Callback(handles.popupmenuAxes,eventdata,handles)
end
set(handles.listboxExistingBuild,'string',char(str))


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

% --- Executes during object creation, after setting all properties.
function lbBuildingClimate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lbBuildingClimate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in lbBuildingVintage.
function lbBuildingVintage_Callback(hObject, eventdata, handles)
% hObject    handle to lbBuildingVintage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

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

global Project COMPONENT Model_dir
Ts = 8760/length(COMPONENT.DemandE);
axes(findobj(gcf,'tag','axes1'))
cla reset
set(gca,'tag','axes1')
hold all

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

val=get(hObject,'value');
str=get(hObject,'string');
currentstr=str{val};
switch currentstr
    case 'Heat Demand'
        PlotType = 'plot';
        Y = COMPONENT.DemandH(1+24/Ts*(day1-1):24/Ts*lastDay);
    case 'Electric Demand'
        PlotType = 'plot';
        Y = COMPONENT.DemandE(1+24/Ts*(day1-1):24/Ts*lastDay);
    case 'Electric Demand by Day'
        PlotType = 'pcolor';
        xgrid=floor(X(1)):floor(X(end));
        ygrid=Ts:Ts:24;
        z=NaN*ones(length(xgrid),length(ygrid));
        for i=1:lastDay-day1+1
            z(i,:)=COMPONENT.DemandE(1+24/Ts*(i-1+day1-1):(i+day1-1)*24/Ts);
        end
    case 'Cooling Demand'
        PlotType = 'plot';
        Y = COMPONENT.DemandC(1+24/Ts*(day1-1):24/Ts*lastDay);
    case 'Electric without HVAC'
        PlotType = 'plot';
        Y = COMPONENT.DemandE(1+24/Ts*(day1-1):24/Ts*lastDay)-COMPONENT.CoolingElectricalLoad(1+24/Ts*(day1-1):24/Ts*lastDay);
end
newlist = {'New Elec Demand';'New Elec w/o HVAC';'New Cooling Demand';};
if max(strcmp(currentstr,newlist))>0
    %load new building profile
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
    buildDir = fullfile(Model_dir,'System Library','Buildings');
    load(strcat(buildDir,filesep,name))
    newBuild = component;
    PlotType = 'plot';
    switch currentstr
        case 'New Elec Demand'
            Y = newBuild.DemandE(1+24/Ts*(day1-1):24/Ts*lastDay);
        case 'New Elec w/o HVAC'
            Y = newBuild.DemandE(1+24/Ts*(day1-1):24/Ts*lastDay)-newBuild.CoolingElectricalLoad(1+24/Ts*(day1-1):24/Ts*lastDay);
        case 'New Cooling Demand'
            Y = newBuild.DemandC(1+24/Ts*(day1-1):24/Ts*lastDay);
    end
end
if strcmp(PlotType,'plot')
    plot(X,Y);
    ylabel('Demand (kW)') 
elseif strcmp(PlotType,'hist')
    hist(Y);
    xlabel('Demand (kW)')
    ylabel('Count')
elseif strcmp(PlotType,'pcolor')
    h=pcolor(xgrid,ygrid,z');
    set(h,'edgecolor','none')
    ylabel('hour')
    colorbar
    ylim([1 24]) 
end
if strcmp(PlotType,'plot') || strcmp(PlotType,'pcolor')
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
    if get(handles.checkboxGeneration,'Value')==1 || (get(handles.checkboxTESshift,'Value')==1 && isfield(Project.System,'TES'))
        ProjectTemp = Project;
        ProjectTemp.Building = COMPONENT;
        Result = runAnalyses(ProjectTemp);
    end
    if get(handles.checkboxGeneration,'Value')==1
        genPow = Result.eOut.GenPower(1+24/Ts*(day1-1):24/Ts*lastDay);
        hold on
        plot(X,genPow,'r');
        hold off
    end
    if get(handles.checkboxTESshift,'Value')==1 && isfield(Project.System,'TES')
        hold on
        plot(X,Result.eOut.DemandEshift(1+24/Ts*(day1-1):24/Ts*lastDay),'g');
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
    steps = 2;
elseif Zoom == 2
    steps = 12;
elseif Zoom == 3
    steps = 52;
elseif Zoom == 4
    steps = 365;
end
set(handles.sliderDate,'Min',1)
set(handles.sliderDate,'Max',steps)
set(handles.sliderDate,'Value',1,'SliderStep',[1/(steps-1),1/(steps-1)])
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
global COMPONENT
DistHeat = COMPONENT.DistHeat;
DistCool = COMPONENT.DistCool;
CHPtemp = COMPONENT.CHPtemp;
NameList = COMPONENT.NameList;
Refrig = COMPONENT.Refrig;
type=get(handles.lbBuildingType, 'value');
climate=get(handles.lbBuildingClimate, 'value');
vintage=get(handles.lbBuildingVintage, 'value');

OpenedRefrigGUI = 0;
vintageName = {'New2010';'New2007';'New2004';'Post1980';'Pre1980';};
climateName = {'_1A_'; '_2A_'; '_2B_'; '_3A_'; '_3B_'; '_3B-Coast_'; '_3C_'; '_4A_'; '_4B_'; '_4C_'; '_5A_'; '_5B_'; '_6A_'; '_6B_'; '_7_'; '_8_';};
buildTypeName = {'SDRest'; 'FFRest'; 'Sch-pri'; 'Sch-sec'; 'LgOff'; 'MdOff'; 'SmOff'; 'MRapt'; 'LgHotel'; 'SmHotel'; 'Hospital'; 'OutP'; 'Retail'; 'StMall'; 'SMarket'; 'ware';};
for i = 1:1:length(type)
    for j = 1:1:length(climate)
        for k = 1:1:length(vintage)
            name = char(strcat(buildTypeName(type(i)),climateName(climate(j)),vintageName(vintage(k))));
            DistHeat(end+1) = get(handles.checkboxDistHeat,'value');
            DistCool(end+1) = get(handles.checkboxDistCool,'value');
            CHPtemp(end+1) = str2double(get(handles.editCHPtemp,'string'));
            NameList{end+1} = cellstr(name); 
            Refrig(end+1).ElecLoad = zeros(8760/4,1);
            Refrig(end).CoolLoad = zeros(8760/4,1);
            if type(i) ==1 || type(i) ==2 || type(i) ==3 || type(i) ==4 || type(i) == 9 || type(i) == 11 || type(i) == 15
                message ={'You have selected a building with considerable refrigeration loads. Would you like to consider heat to cooling applications?'};
                choice = menu(message,'Yes','No');
                if choice ==1
                    OpenedRefrigGUI = 1;
                    Refrig(end) = refrigerationChooser(name,DistCool);
                    uiwait(gcf)
                end
            end
        end
    end
end
reloadBuildings(DistHeat,DistCool,CHPtemp,NameList,Refrig)
listboxExistingBuild_Callback(handles.listboxExistingBuild, eventdata, handles)
if OpenedRefrigGUI
    pushbuttonDone_Callback(hObject, eventdata, handles)
end
% --- Executes on button press in pushbuttonRemove.
function pushbuttonRemove_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonRemove (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global COMPONENT
val = get(handles.listboxExistingBuild,'value');
DistHeat = [];
DistCool = [];
CHPtemp = [];
NameList = {};
Refrig={};
n = length(COMPONENT.NameList);
for i = 1:1:n
    if i ~=val
        NameList{end+1} = COMPONENT.NameList{i};
        DistHeat(end+1) = COMPONENT.DistHeat(i);
        DistCool(end+1) = COMPONENT.DistCool(i);
        CHPtemp(end+1) = COMPONENT.CHPtemp(i);
        if i==1 || (i==2&&val==1)
            Refrig = COMPONENT.Refrig(i);
        else Refrig(end+1) = COMPONENT.Refrig(i);
        end
    end
end
reloadBuildings(DistHeat,DistCool,CHPtemp,NameList,Refrig)
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
COMPONENT.CHPtemp = [];
COMPONENT.Refrig ={};
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

% --- Executes on button press in pushbuttonResults.
function pushbuttonResults_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonResults (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close(gcf)
ViewResults2(1)

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


function editCHPtemp_Callback(hObject, eventdata, handles)
% hObject    handle to editCHPtemp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function editCHPtemp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editCHPtemp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function reloadBuildings(DistHeat,DistCool,CHPtemp,NameList,Refrig)
global COMPONENT Model_dir
names = fieldnames(COMPONENT);
COMPONENT = rmfield(COMPONENT,names);
%%double -check , rebuild demand profiles from scratch
buildDir = fullfile(Model_dir,'System Library','Building');
COMPONENT.DemandE = zeros(35040,1);
COMPONENT.DemandC = zeros(35040,1);
COMPONENT.DemandH = zeros(35040,1);
COMPONENT.NonDistH = zeros(35040,1);
COMPONENT.CoolingElectricalLoad = zeros(35040,1);
COMPONENT.Refrig = Refrig;
COMPONENT.DistCool = DistCool;
COMPONENT.DistHeat = DistHeat;
COMPONENT.CHPtemp = CHPtemp;
COMPONENT.NameList = NameList;
if ~isfield(COMPONENT,'NameList')
    name = char(COMPONENT.Name);
    load(fullfile(buildDir,'RealBuildingData',name))
    COMPONENT.DemandE = component.DemandE;
    COMPONENT.DemandC = component.DemandC;
    COMPONENT.DemandH = component.DemandH;
    COMPONENT.CoolingElectricalLoad = DistCool*component.CoolingElectricalLoad;
    COMPONENT.NonDistH = (1-DistHeat)*component.DemandH;
else
    for i = 1:1:length(COMPONENT.NameList)
        name = char(COMPONENT.NameList{i});
        load(fullfile(buildDir,name))
        %%% edit exterior lighting profile
            AbsMin =min(component.DemandE);
            LightLoad = max(component.ExteriorLight);
            AddLight = (component.ExteriorLight==0).*(component.DemandE<=(AbsMin+LightLoad));
            component.DemandE = component.DemandE+AddLight*LightLoad;
%             COMPONENT.ExteriorLight = COMPONENT.ExteriorLight+ component.ExteriorLight+AddLight*LightLoad;
        %%% 
        if isfield(Refrig(i),'ElecLoad') && length(Refrig(i).ElecLoad)==length(component.DemandE) && max(Refrig(i).ElecLoad)>0
            COMPONENT.DemandE = COMPONENT.DemandE + component.DemandE - Refrig(i).ElecLoad;
            COMPONENT.DemandC = COMPONENT.DemandC + DistCool(i)*component.DemandC + Refrig(i).CoolLoad;
        else COMPONENT.DemandE = COMPONENT.DemandE + component.DemandE;
            COMPONENT.DemandC = COMPONENT.DemandC + DistCool(i)*component.DemandC;
        end
        COMPONENT.DemandH = COMPONENT.DemandH + DistHeat(i)*component.DemandH;
        COMPONENT.NonDistH = COMPONENT.NonDistH + (1-DistHeat(i))*component.DemandH;
        COMPONENT.CoolingElectricalLoad = COMPONENT.CoolingElectricalLoad + DistCool(i)*component.CoolingElectricalLoad;
    end
end

