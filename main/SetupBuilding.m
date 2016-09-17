function SetupBuilding(varargin)
% SETUPBUILDING MATLAB code for SetupBuilding.fig
%      SETUPBUILDING, by itself, creates a new SETUPBUILDING or raises the existing
%      singleton*.
%
%      H = SETUPBUILDING returns the handle to a new SETUPBUILDING or the handle to
%      the existing singleton*.
%
%      SETUPBUILDING('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SETUPBUILDING.M with the given input arguments.
%
%      SETUPBUILDING('Property','Value',...) creates a new SETUPBUILDING or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SetupBuilding_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SetupBuilding_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SetupBuilding

% Last Modified by GUIDE v2.5 09-Apr-2015 09:34:50

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SetupBuilding_OpeningFcn, ...
                   'gui_OutputFcn',  @SetupBuilding_OutputFcn, ...
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


% --- Executes just before SetupBuilding is made visible.
function SetupBuilding_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SetupBuilding (see VARARGIN)

% Choose default command line output for SetupBuilding
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
global update Plant figHandle BUILDING DataOriginal Model_dir
DataOriginal=Plant.Data;
figHandle = handles.figure1;
update = 0;
if isfield(Plant.Data,'NameList')
    BUILDING.NameList = Plant.Data.NameList;
    if isfield(Plant.Data.Demand,'E')
        BUILDING.DemandE = Plant.Data.Demand.E;
    end
    if isfield(Plant.Data.Demand,'H')
        BUILDING.DistHeat(1) = 1;
        BUILDING.DemandH = Plant.Data.Demand.H;
    else BUILDING.DistHeat(1) = 0;
    end
    if isfield(Plant.Data.Demand,'C')
        BUILDING.DistCool(1) = 1;
        BUILDING.DemandC = Plant.Data.Demand.C;
    else BUILDING.DistCool(1) = 0;
    end
else
    BUILDING = [];
end
BuildType = {'Restaurant: full-service (sit down)';'Restaurant: quick-service (fast food)';'School: primary school';'School: secondary school';'Office: large office';'Office: medium office';'Office: small office';'Mid-rise apartment building';'Hospitality: large hotel';'Hospitality: small hotel/motel';'Health care: large hospital';'Health care: outpatient facility';'Retail: big-box, standalone retail store';'Retail: retail store located in a strip mall';'Retail: supermarket';'Unrefrigerated warehouse';};
Climate = {'Miami (ASHRAE 1A)';'Houston (ASHRAE 2A)';'Phoenix (ASHRAE 2B)';'Atlanta (ASHRAE 3A)';'Las Vegas (ASHRAE 3B-Inland)';'Los Angeles (ASHRAE 3B-Coast)';'San Francisco (ASHRAE 3C)';'Baltimore (ASHRAE 4A)';'Albuquerque (ASHRAE 4B)';'Seattle (ASHRAE 4C)';'Chicago (ASHRAE 5A)';'Boulder (ASHRAE 5B)';'Minneapolis (ASHRAE 6A)';'Helena, MT (ASHRAE 6B)';'Duluth, MN (ASHRAE 7)';'Fairbanks, AK (ASHRAE 8)';};
Vintage = {'2010 construction (ASHRAE 90.1-2010)';'2007 construction (ASHRAE 90.1-2007)';'2004 construction 90.1-2004';'“Post-1980” construction (ASHRAE 90.1-1989)';'“Pre-1980” construction';};
%Profiles = {'Demand C'; 'Demand C2'; 'Demand E'; 'Demand E C'; 'Demand E H'; 'Demand E2'; 'Demand E3';};
set(handles.lbBuildingType, 'string', BuildType(:,1), 'Max', length(BuildType(:,1)), 'Min', 1)
set(handles.lbBuildingClimate, 'string', Climate(:,1), 'Max', length(Climate(:,1)), 'Min', 1)
set(handles.lbBuildingVintage, 'string', Vintage(:,1), 'Max', length(Vintage(:,1)), 'Min', 1)
%set(handles.lbBuildingProfiles, 'string', Profiles(:,1), 'Max', length(Profiles(:,1)), 'Min', 1)

MATfiles=dir(fullfile(Model_dir,'component library','Weather','*.mat'));
Climate ={};
for i = 1:1:length(MATfiles)
    Climate(end+1) = cellstr(strrep(MATfiles(i).name,'.mat',''));
end
%Plant data is an empty vector until it is specified by user
Plant.Data.Weather = zeros(1,length(Climate));
val = find(strcmp(Climate,Plant.Data.Weather));
if isempty(val)
    val=1;
end
set(handles.popupmenuTempProfile, 'string', Climate,'Value', val)
listboxExistingBuild_Callback(handles.listboxExistingBuild,eventdata,handles)

set(handles.checkboxDistHeat,'value',0);    
set(handles.checkboxDistCool,'value',0); 
days = round(Plant.Data.Timestamp(end) - Plant.Data.Timestamp(1));
set(handles.sliderZoom,'Min',1,'Max',4,'Value',1,'SliderStep',[1/3,1/3])
set(handles.sliderDate,'Min',1,'Max',2,'Value',1,'SliderStep',[1/(days-1),1/(days-1)])
set(handles.sliderDate,'Max',2)
sliderZoom_Callback(handles.sliderZoom, eventdata, handles)

% MATfilesP=dirfullfile(Model_dir, 'Plant', 'Load Profiles', '*.mat');
% Profiles ={};
% for i = 1:1:length(MATfilesP)
%     Profiles(end+1) = cellstr(strrep(MATfiles(i).name, '.mat', ''));
% end

if isfield(Plant.Data,'HistProf')&& ~isempty(Plant.Data.HistProf)
    list = {};
    if isfield(Plant.Data.HistProf,'Temperature')
        list(end+1) = {'Temperature'};
    end
    if isfield(Plant.Data.HistProf,'E')
        list(end+1) = {'Electric'};
    end
    if isfield(Plant.Data.HistProf,'H')
        list(end+1) = {'Heating'};
    end
    if isfield(Plant.Data.HistProf,'C')
        list(end+1) = {'Cooling'};
    end
    set(handles.popupmenuHistFits,'string',list,'value',1);
    plotHistFits(handles,2)
else
    set(handles.popupmenuHistFits,'string',{'Temperature';},'value',1);    
end


% --- Outputs from this function are returned to the command line.
function varargout = SetupBuilding_OutputFcn(hObject, eventdata, handles) 
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
global BUILDING Plant
vintageName = {'New2010';'New2007';'New2004';'Post1980';'Pre1980';};
climateName = {'_1A_'; '_2A_'; '_2B_'; '_3A_'; '_3B_'; '_3B-Coast_'; '_3C_'; '_4A_'; '_4B_'; '_4C_'; '_5A_'; '_5B_'; '_6A_'; '_6B_'; '_7_'; '_8_';};
buildTypeName = {'SDRest'; 'FFRest'; 'Sch-pri'; 'Sch-sec'; 'LgOff'; 'MdOff'; 'SmOff'; 'MRapt'; 'LgHotel'; 'SmHotel'; 'Hospital'; 'OutP'; 'Retail'; 'StMall'; 'SMarket'; 'ware';};
profiles = {'Demand_C'; 'Demand_C2'; 'Demand_E'; 'Demand_E_C'; 'Demand_E_H'; 'Demand_E2'; 'Demand_E3';};

BuildType = get(handles.lbBuildingType, 'string');
Climate = get(handles.lbBuildingClimate, 'string');
Vintage = get(handles.lbBuildingVintage, 'string');
plotlist={'Temperature'};
    %build list of existing names
str={};
if ~isempty(BUILDING.NameList) && ~isfield(BUILDING,'DemandE')
    if isfield(Plant.Data.Demand,'E')
        BUILDING.DemandE = Plant.Data.Demand.E;
    end
    if isfield(Plant.Data.Demand,'H')
        BUILDING.DistHeat(1) = 1;
        BUILDING.DemandH = Plant.Data.Demand.H;
    else BUILDING.DistHeat(1) = 0;
    end
    if isfield(Plant.Data.Demand,'C')
        BUILDING.DistCool(1) = 1;
        BUILDING.DemandC = Plant.Data.Demand.C;
    else BUILDING.DistCool(1) = 0;
    end

elseif isempty(BUILDING.NameList)
    str(end+1) = cellstr('None Loaded');
    BUILDING.DistHeat = [];
    BUILDING.DistCool = [];
%     BUILDING.DemandE = [];%zeros(96*365,1);
%     BUILDING.DemandC = [];%BUILDING.DemandE;
%     BUILDING.DemandH = [];%BUILDING.DemandE;
end
for i = 1:1:length(BUILDING.NameList)
    name = char(BUILDING.NameList{i});
    ind =strfind(name,'_'); %the index of the _ in the name
    if isempty(ind) || isempty(find(strcmp(name(1:ind(1)-1),buildTypeName)))
        str(end+1) = BUILDING.NameList(i);
        BUILDING.DistHeat(i) = BUILDING.DistHeat(1);
        BUILDING.DistCool(i) = BUILDING.DistCool(1);
    else
        build = find(strcmp(name(1:ind(1)-1),buildTypeName));
        clim = find(strcmp(name(ind(1):ind(2)),climateName));
        %vin = find(strcmp(name(ind(2)+1:ind(3)-1),vintageName)); use this
        %once district heating options are implemented
        vin = find(strcmp(name(ind(2)+1:length(name)),vintageName));
        %uncomment the below once the district heating and cooling options
        %are implemented
%         if strcmp(name(ind(3)+1),'h')
%             BUILDING.DistHeat(i) = 1;
%         else BUILDING.DistHeat(i) = 0;
%         end
%         if strcmp(name(ind(4)+1),'c')
%             BUILDING.DistCool(i) = 1;
%         else BUILDING.DistCool(i) = 0;
%         end
        str(end+1) = {char(strcat(BuildType(build), Climate(clim), Vintage(vin)))};
    end
end
if isfield(BUILDING,'DemandE')
    plotlist(end+1) = {'Electric Demand'};
end
if isfield(BUILDING,'DemandC')
    plotlist(end+1) = {'Cooling Demand'};
end
if isfield(BUILDING,'DemandH')
    plotlist(end+1) = {'Heating Demand'};
end

set(handles.popupmenuAxes,'string',plotlist,'value',1)
set(handles.listboxExistingBuild,'string',char(str))


% --- Executes during object creation, after setting all properties.
function listboxExistingBuild_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --- Executes on selection change in lbBuildingType.
function lbBuildingType_Callback(hObject, eventdata, handles)
popupmenuAxes_Callback(hObject, eventdata, handles)
% --- Executes during object creation, after setting all properties.
function lbBuildingType_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --- Executes on selection change in lbBuildingClimate.
function lbBuildingClimate_Callback(hObject, eventdata, handles)
% --- Executes during object creation, after setting all properties.
popupmenuAxes_Callback(hObject, eventdata, handles)
function lbBuildingClimate_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --- Executes on selection change in lbBuildingVintage.
function lbBuildingVintage_Callback(hObject, eventdata, handles)
% --- Executes during object creation, after setting all properties.
popupmenuAxes_Callback(hObject, eventdata, handles)
function lbBuildingVintage_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pushbuttonDone.
function pushbuttonDone_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonDone (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Plant BUILDING update figHandle
Plant.Data.NameList = BUILDING.NameList;
if isfield(BUILDING,'DemandE')
    if isfield(BUILDING,'DemandE')
        Plant.Data.Demand.E = BUILDING.DemandE;
    end
    if isfield(BUILDING,'DemandC')
        Plant.Data.Demand.C = BUILDING.DemandC;
    end
    if isfield(BUILDING,'DemandH')
        Plant.Data.Demand.H = BUILDING.DemandH;
    end
end
if update ==0
    h = menu('Historical Demand Curve Fits Not Updated ','Update Now','Do not update');
    if h ==1
        pushbuttonHistFit_Callback(hObject, eventdata, handles,1)
    end
end
close(figHandle)


function popupmenuAxes_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuAxes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Plant BUILDING Model_dir
val=get(handles.popupmenuAxes,'value');
list=get(handles.popupmenuAxes,'string');
currentstr = char(list(val));
PlotIndex = PlotWindow2(handles, Plant.Data.Timestamp);
X = Plant.Data.Timestamp(PlotIndex(1):PlotIndex(2))';
axes(handles.axes1)
if strcmp(char(list(val)),'Temperature')
    A = Plant.Data.Temperature;
    Y = A(PlotIndex(1):PlotIndex(2))';
    if isstruct(Plant.Data.HistProf) %if you have already loaded the historical profiles, the structure of HistProf is different
        PlotData(X,Y,0,'T',Plant.Data,[],[], [])
    else
        PlotData(X,Y,0,'T',Plant.Data,Plant.Data.HistProf,[], [])
    end
else
    A = BUILDING.(strcat('Demand',currentstr(1)));
    Y = A(PlotIndex(1):PlotIndex(2))';
    
    %load from selected building
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
    load(fullfile(Model_dir, 'component library','Buildings',name));
    A = component.(strcat('Demand',currentstr(1)));
%     D = datevec(PlotIndex(2));%PlotIndex(2) gives Plant.Data.Timestamp
%     which is the number of steps per day, not the date to start on
    D = datevec(X(1));
    D2 = datenum([D(1) 1 1 0 0 0])+linspace(0,365,96*365);%from 0 to end of year with steps of 15 min
   %take all the indicies in A that align with the timestamps in X. X is 96
   %long (one day), but A is 35040 (whole year)
    Xi = nnz(D2<=X(1));%find the index of A that corresponds with start time X
    Xend = nnz(D2<=X(end));%find the index of A that corresponds with the end time X
    Y(2,:) = A(Xi:Xend)';%use the data for that day
%     for t = 1:1:length(X)
%         Xi = nnz(D2<=X(t));
%         Y(t,2) = A(Xi);
%     end
    %plot    
    if isstruct(Plant.Data.HistProf) %if you have loaded the historical data, then the structure of Plant.Data.HistProf is different
        PlotData(X,Y,0,currentstr(1),Plant.Data,[],[],[])
    else PlotData(X,Y,0,currentstr(1),Plant.Data,Plant.Data.HistProf,[],[])
    end
    legend('Current Microgrid','Selected Building')
end


function [PlotIndex] = PlotWindow2(handles, timestamp)
a = datevec(timestamp(1));
b = datevec(timestamp(end));
D1 = datenum([a(1) a(2) a(3)]);
months = max(1,12*(b(1)-a(1))+b(2)-a(2));
years = max(1,b(1)-a(1));
days = round(timestamp(end)-timestamp(1));
weeks = ceil(days/7);
Zoom = max(1,round(get(handles.sliderZoom,'Value')));
if Zoom == 1
    day = 1+round((get(handles.sliderDate,'Value')-1)*(days-1));
    endday = day;
elseif Zoom ==2
    day = 1+7*round((get(handles.sliderDate,'Value')-1)*(weeks-1));
    endday = day+6;
elseif Zoom ==3
    month = 1+ round((get(handles.sliderDate,'Value')-1)*(months-1));
    day = datenum(a(1),month,1)-datenum(a(1),1,1)+1;
    endday = day+(datenum(a(1),month+1,1)-datenum(a(1),month,1)-1);
elseif Zoom ==4
    year = a(1)+ round((get(handles.sliderDate,'Value')-1)*(years-1));
    day = 1+round(datenum(year,1,1)-round(timestamp(1)));
    endday = round(datenum(year+1,1,1)-round(timestamp(1)));
end
PlotIndex(1) = max(1,nnz(timestamp<=(D1+(day-1))));
PlotIndex(2) = nnz(timestamp<=(D1+endday));
if PlotIndex(2)>length(timestamp)
    set(handles.sliderZoom,'Value',Zoom-1)
    PlotWindow2(handles, timestamp)
end
    

% --- Executes on slider movement.
function sliderZoom_Callback(hObject, eventdata, handles)
% hObject    handle to sliderZoom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
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

% --- Executes on button press in pushbuttonAdd.
function pushbuttonAdd_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonAdd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global BUILDING Model_dir update
update =0;
type=get(handles.lbBuildingType, 'value');
climate=get(handles.lbBuildingClimate, 'value');
vintage=get(handles.lbBuildingVintage, 'value');

vintageName = {'New2010';'New2007';'New2004';'Post1980';'Pre1980';};
climateName = {'_1A_'; '_2A_'; '_2B_'; '_3A_'; '_3B_'; '_3B-Coast_'; '_3C_'; '_4A_'; '_4B_'; '_4C_'; '_5A_'; '_5B_'; '_6A_'; '_6B_'; '_7_'; '_8_';};
buildTypeName = {'SDRest'; 'FFRest'; 'Sch-pri'; 'Sch-sec'; 'LgOff'; 'MdOff'; 'SmOff'; 'MRapt'; 'LgHotel'; 'SmHotel'; 'Hospital'; 'OutP'; 'Retail'; 'StMall'; 'SMarket'; 'ware';};
for i = 1:1:length(type)
    for j = 1:1:length(climate)
        for k = 1:1:length(vintage)
            BUILDING.DistHeat(end+1) = get(handles.checkboxDistHeat,'value');
            BUILDING.DistCool(end+1) = get(handles.checkboxDistCool,'value');
            if BUILDING.DistHeat(end)==1
                h='h';
            else h = 'o';
            end
            if BUILDING.DistCool(end)==1
                c='c';
            else c = 'o';
            end
            name = char(strcat(buildTypeName(type(i)),climateName(climate(j)),vintageName(vintage(k))));% haven't incorporated different buildings for heating and cooling options yet'_',h,'_',c));
            BUILDING.NameList{end+1} = cellstr(name); 
            if isfield(BUILDING,'DemandE') 
                E = BUILDING.DemandE;
            else E =0;
            end
            if isfield(BUILDING,'DemandC') 
                C = BUILDING.DemandC;
            else C =0;
            end
            if isfield(BUILDING,'DemandH') 
                H = BUILDING.DemandH;
            else H =0;
            end
            load(fullfile(Model_dir, 'component library','Buildings',name));
            BUILDING.DemandE = E + component.DemandE - BUILDING.DistCool(end)*component.CoolingElectricalLoad;
            if BUILDING.DistCool(end)==1
                BUILDING.DemandC = C + BUILDING.DistCool(end)*component.DemandC;
%             else
%                 BUILDING.DemandC = C + component.DemandC;
            end
            if BUILDING.DistHeat(end)==1
                BUILDING.DemandH = H + BUILDING.DistHeat(end)*component.DemandH;
            %else
               % BUILDING.DemandH = H + component.DemandH;
            end
        end
    end
end
listboxExistingBuild_Callback(handles.listboxExistingBuild, eventdata, handles)
popupmenuAxes_Callback(hObject, eventdata, handles)

% --- Executes on button press in pushbuttonRemove.
function pushbuttonRemove_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonRemove (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global BUILDING Model_dir update
update =0;
val = get(handles.listboxExistingBuild,'value');
buildTypeName = {'SDRest'; 'FFRest'; 'Sch-pri'; 'Sch-sec'; 'LgOff'; 'MdOff'; 'SmOff'; 'MRapt'; 'LgHotel'; 'SmHotel'; 'Hospital'; 'OutP'; 'Retail'; 'StMall'; 'SMarket'; 'ware';};
DistHeat = [];
DistCool = [];
NameList = {};
n = length(BUILDING.NameList);
name = char(BUILDING.NameList{val});
ind =strfind(name,'_');
build = find(strcmp(name(1:ind(1)-1),buildTypeName),1);

if isempty(build) %not a building from library
%     BUILDING.DemandE = 0*BUILDING.DemandE;
%     BUILDING.DemandC = BUILDING.DemandE;
%     BUILDING.DemandH = BUILDING.DemandE;
    for i = 1:1:n %start from scratch and add all other buildings from library
       if i ~=val
            name = char(BUILDING.NameList{i});
            build = find(strcmp(name(1:ind(1)-1),buildTypeName),1);
            if ~isempty(build)
                NameList(end+1) = BUILDING.NameList{i};
                DistHeat(end+1) = BUILDING.DistHeat(i);
                DistCool(end+1) = BUILDING.DistCool(i);
                if isfield(BUILDING,'DemandE') 
                    E = BUILDING.DemandE;
                else E =0;
                end
                if isfield(BUILDING,'DemandC') 
                    C = BUILDING.DemandC;
                else C =0;
                end
                if isfield(BUILDING,'DemandH') 
                    H = BUILDING.DemandH;
                else H =0;
                end
                load(fullfile(Model_dir, 'component library','Buildings',char(NameList(i))));
                BUILDING.DemandE = E + component.DemandE - BUILDING.DistCool(end)*component.CoolingElectricalLoad;
                if BUILDING.DistCool(end)==1
                    BUILDING.DemandC = C + BUILDING.DistCool(end)*component.DemandC;
                end
                if BUILDING.DistHeat(end)==1
                    BUILDING.DemandH =H + BUILDING.DistHeat(end)*component.DemandH;
                end
            end
       end
    end
else %Remove the building that was in the library
    if isfield(BUILDING,'DemandE') 
        E = BUILDING.DemandE;
    else E =0;
    end
    if isfield(BUILDING,'DemandC') 
        C = BUILDING.DemandC;
    else C =0;
    end
    if isfield(BUILDING,'DemandH') 
        H = BUILDING.DemandH;
    else H =0;
    end
    load(fullfile(Model_dir, 'component library','Buildings',name));
    BUILDING.DemandE = E -(component.DemandE - BUILDING.DistCool(end)*component.CoolingElectricalLoad);
    if BUILDING.DistCool(val)==1
        BUILDING.DemandC = C - BUILDING.DistCool(val)*component.DemandC;
    end
    if BUILDING.DistHeat(val)==1
        BUILDING.DemandH = H - BUILDING.DistHeat(val)*component.DemandH;
    end
    for i = 1:1:n
        if i ~=val
            NameList{end+1} = BUILDING.NameList{i};
            DistHeat(end+1) = BUILDING.DistHeat(i);
            DistCool(end+1) = BUILDING.DistCool(i);
        end
    end
end
BUILDING.DistHeat = DistHeat;
BUILDING.DistCool = DistCool;
BUILDING.NameList = NameList;
set(handles.listboxExistingBuild,'Value',1)
listboxExistingBuild_Callback(handles.listboxExistingBuild, eventdata, handles)
popupmenuAxes_Callback(hObject, eventdata, handles)

% --- Executes on button press in pushbuttonClear.
function pushbuttonClear_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonClear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global BUILDING update
update =0;
clear BUILDING
BUILDING.NameList = {};
BUILDING.DistHeat = [];
BUILDING.DistCool = [];
% BUILDING.DemandE = 0*BUILDING.DemandE;
% BUILDING.DemandC = BUILDING.DemandE;
% BUILDING.DemandH = BUILDING.DemandE;
set(handles.listboxExistingBuild,'Value',1)
listboxExistingBuild_Callback(handles.listboxExistingBuild, eventdata, handles)
popupmenuAxes_Callback(hObject, eventdata, handles)

% --- Executes on button press in pushbuttonSwitch.
function pushbuttonSwitch_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonSwitch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
pushbuttonRemove_Callback(hObject, eventdata, handles)
pushbuttonAdd_Callback(hObject, eventdata, handles)


% --- Executes on button press in checkboxDistHeat.
function checkboxDistHeat_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxDistHeat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global BUILDING
BUILDING.DistHeat = 1;

% --- Executes on button press in checkboxDistCool.
function checkboxDistCool_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxDistCool (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global BUILDING
BUILDING.DistCool = 1;


% --- Executes on button press in pushbuttonHistFit.
function pushbuttonHistFit_Callback(hObject, eventdata, handles,skipPlot)
% hObject    handle to pushbuttonHistFit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global update Plant BUILDING
update =1;
list = {'Temperature'};
if isfield(BUILDING,'DemandE')
    Plant.Data.Demand.E = BUILDING.DemandE;
    list(end+1) = {'Electricity'};
elseif isfield(Plant.Data.Demand,'E')
    rmfield(Plant.Data.Demand,'E')
end
if isfield(BUILDING,'DemandC')
    Plant.Data.Demand.C = BUILDING.DemandC;
    list(end+1) = {'Cooling'};
elseif isfield(Plant.Data.Demand,'C')
    rmfield(Plant.Data.Demand,'C')
end
if isfield(BUILDING,'DemandH')
    Plant.Data.Demand.H = BUILDING.DemandH;
    list(end+1) = {'Heating'};
elseif isfield(Plant.Data.Demand,'H')
    rmfield(Plant.Data.Demand,'H')
end
list(end+1) = cellstr('All');
[s,v] = listdlg('PromptString','Choose Variable to Fit', 'SelectionMode','single','ListString',list);
if v ==1 
    name = char(list(s));
    if strcmp(name,'All')
        s = 1;
        calculateFit('T')
        n = char(fieldnames(Plant.Data.Demand));
        for i = 1:1:length(n)
            if max(Plant.Data.Demand.(n(i)))>0
                calculateFit(n(i))
            end
        end
    else calculateFit(name(1))
    end
end

set(handles.popupmenuHistFits,'string',list(1:end-1),'value',s);
if ~exist('skipPlot','var') || isempty(skipPlot)
    plotHistFits(handles,0)
else plotHistFits(handles,skipPlot)
end

% --- Executes on selection change in popupmenuHistFits.
function popupmenuHistFits_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuHistFits (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
plotHistFits(handles,0)

% --- Executes during object creation, after setting all properties.
function popupmenuHistFits_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuHistFits (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function sliderZoom2_Callback(hObject, eventdata, handles)
% hObject    handle to sliderZoom2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
demand = get(handles.popupmenuHistFits, 'Value');
if demand == 1%if plotting temp, there is no ZLim, use y
    a = get(handles.axes2, 'YLim');
elseif demand == 2 %if plotting electric, use ZLim from 0 to top
    a = [0, mean(get(handles.axes2, 'ZLim'))*2];
elseif demand == 3 %if plotting heat, use Zlim
    a = [-mean(get(handles.axes2,'ZLim'))*2, mean(get(handles.axes2, 'ZLim'))*4];%typically from -50 to 100, or -100 to 200 etc.
end
Zoom = (1-get(handles.sliderZoom2,'Value'));
if Zoom == 0 %the axis can never show a range of zero
    Zoom = 0.001;
end
b(1) = mean(a)-Zoom*(mean(a)-a(1));
b(2) = mean(a)+Zoom*(a(2)-mean(a));
if demand ==1
    set(handles.axes2,'YLim',b)
else set(handles.axes2, 'ZLim',b)
end

% --- Executes during object creation, after setting all properties.
function sliderZoom2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderZoom2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


function calculateFit(name)
global Plant HistProf
h=waitbar(0,'Recalculating surface fit');
Steps = round(1/(Plant.Data.Timestamp(2)-Plant.Data.Timestamp(1)));%points in 1 day of data
Resolution =24/Steps;
[~,m,~,hour,minutes] = datevec(Plant.Data.Timestamp(2:end-1));%avoid checking month 1st and last point if starting or ending at 00
m = [m(1) m m(end)];%make m the same lenght as data
hour = [max(0,floor(hour(1)-Resolution)) hour min(24,hour(end)+floor((minutes(end)+60*Resolution)/60))]';%make hour the same lenght as data
minutes = [max(0,minutes(1)-60*Resolution) minutes mod(minutes(end)+60*Resolution,60)]';%make hour the same lenght as data
months = unique(m);
Xs = 1+round(Steps*(ceil(Plant.Data.Timestamp(1))-Plant.Data.Timestamp(1)));
if strcmp(name,'T') 
    Z = ones(12,24);
    for i = 1:1:length(months)
        waitbar(i/12,h,'Recalculating surface fit');
        D = datevec(Plant.Data.Timestamp(Xs));
        days = floor(nnz(m==D(2))/Steps); %data points in month/points per day
        Total = zeros(24,1);
        points = zeros(24,1);
        for d = 1:1:days
            for j = 1:1:24
                XFs = Xs+round(Steps/24)-1;
                Total(j) = Total(j)+sum(Plant.Data.Temperature(Xs:XFs));
                points(j) = points(j) +(XFs-Xs+1);
                Xs = XFs+1;
            end
        end
        Z(D(2),:) = (Total./points)';
    end
    if length(months)<12
        if nnz(months==1)==1 && nnz(months==12)==1
            j =1;
            for i= 2:1:11
                if nnz(months==i)==1
                    j=j+1;
                else Z(i,:) = Z(j,:);
                end
            end
        else
            for i= 1:1:min(months)-1
                Z(i,:) = Z(min(months),:);
            end
            for i= max(months)+1:1:12
                Z(i,:) = Z(max(months),:);
            end
        end
    end
    Plant.Data.HistProf.Temperature = Z;
else %this should make a surface fit for Demand.E, Demand.H, and Demand.C
    monthNames = cellstr({'Jan';'Feb';'Mar';'Apr';'May';'Jun';'Jul';'Aug';'Sep';'Oct';'Nov';'Dec';});
    dateDay = round(Plant.Data.Timestamp(2:end-1));
    for i = 1:1:length(dateDay)
    BusDay(i) = isbusday(dateDay(i),Plant.Data.Holidays,[1 0 0 0 0 0 1]);
    end
    BusDay = [BusDay(1); transpose(BusDay); BusDay(end)];%make m the same lenght as data
    if length(months)<=4 %make 1 surface fit
        Z = Plant.Data.Demand.(name(1))';
        Y = Plant.Data.Temperature;
        X = hour+minutes/60;
        n =length(X);
        if Steps>100 % too many data points for a surface fit, average each hour
            r = ceil(Steps/24);
            X2 = zeros(3000,1);
            Y2 = zeros(3000,1);
            Z2 = zeros(3000,1);
            for i=1:1:floor(n/r)%turn while loop into for loop to make loop faster
%             i=1;
%             while i*r<n
                X2(i)=sum(X((i-1)*r+1:i*r))/r;
                Y2(i)=sum(Y((i-1)*r+1:i*r))/r;
                Z2(i)=sum(Z((i-1)*r+1:i*r))/r;
%                 i=i+1;
            end
            X = X2(1:i-1);%remove indicies that aren't used in the loop
            Y = Y2(1:i-1);
            Z = Z2(1:i-1);
        end
        Plant.Data.HistProf.(name(1)) = fit([X Y],Z,'lowess');
    else %split surface fit by month and weekday/weekend
        for i = 1:1:length(months)
            waitbar((i-1)/length(months),h,strcat('Recalculating surface fit for ',name,' for the month of ',char(monthNames(months(i)))));
            Z = Plant.Data.Demand.(name(1));
            Y = Plant.Data.Temperature;
            X = hour+minutes/60;
%             waitbar(i/length(months),h,strcat('Recalculating surface fit for ',name,'for the month of ',char(monthNames(months(i)))));
%             Z = Plant.Data.Demand.(name(1)).*(m==months(i))';
%             Y = Plant.Data.Temperature.*(m==months(i))';
%             X = hour.*(m==months(i))'+minutes/60.*(m==months(i))';
            if Steps>150 % too many data points for a surface fit, average each 10 min
                r = ceil(Steps/(24*6));
                n =length(X);
%                 j=1;
                X2=zeros(2544);
                Y2=zeros(2544);
                Z2=zeros(2544);
                for j=1:1:floor(n/r)%turn the while loop into a for loop to make loop faster
%                 while j*r<n
                    X2(j)=sum(X((j-1)*r+1:j*r))/r;
                    Y2(j)=sum(Y((j-1)*r+1:j*r))/r;
                    Z2(j)=sum(Z((j-1)*r+1:j*r))/r;
%                     j=j+1;
                end
                X = X2(1:j-1);%remove the indicies that aren't used.
                Y = Y2(1:j-1);
                Z = Z2(1:j-1);
            end
            % business days
            H1 = BusDay((find(m==months(i))));%find which days are businessdays
            H1non0=find(H1);%find the indexes of all the nonzeros (all the businessdays are nonzero)
            Z1 = Z(H1non0);%remove the days that are non business days
            Y1 = Y(H1non0);
            X1 = X(H1non0);%cannot use nonzeros function because there are some elements in X1 that are 0
            % weekends/holidays
            H2 = 1-H1;
            H2non0=find(H2);
            Z2 = Z(H2non0);
            Y2 = Y(H2non0);
            X2 = X(H2non0);
            Plant.Data.HistProf.(name(1)).(strcat(char(monthNames(months(i))),'WeekDay')) = fit([X1, Y1],Z1,'lowess');%'poly23');%'loess');
            Plant.Data.HistProf.(name(1)).(strcat(char(monthNames(months(i))),'WeekEnd')) = fit([X2, Y2],Z2,'lowess');%'poly23');%'loess');
        end
    end
end
HistProf = Plant.Data.HistProf;
close(h)

function plotHistFits(handles,skip)
if skip==0 || skip ==2
    global HistProf
    hold off
    list = get(handles.popupmenuHistFits,'string');
    val=get(handles.popupmenuHistFits,'value');
    sel = char(list{val});

    if strcmp(sel,'Temperature')
        A = HistProf.Temperature;
        B = {'All';'Jan';'Feb';'Mar';'Apr';'May';'Jun';'Jul';'Aug';'Sep';'Oct';'Nov';'Dec';};
        if skip ==2
            s=1;
        else
            [s,v] = listdlg('PromptString','Choose a profile to plot', 'SelectionMode','single','ListString',B);
        end
        if s==1
            Y = [A(:,end),A];
        else Y = [A(s-1,end),A(s-1,:)];
        end
        X = linspace(0,24,length(Y(1,:)));
        axes(handles.axes2)
        hold off
        plot(X',Y')
        ylabel('Temperature')
    else
    A = HistProf.(sel(1));
        B = fields(A);
        if isempty(B)
            plot(A)
            ylabel(char(sel));
        else
         [s,v] = listdlg('PromptString','Choose a profile to plot', 'SelectionMode','single','ListString',B);
            h=waitbar(.5,'Plotting Historical Fit');
            axes(handles.axes2)
            hold off
            plot(A.(char(B(s))))
            ylabel(char(sel))
            close(h)
        end
    end
    xlim([0,24])
    set(handles.axes2,'XTick',linspace(0,24,7))
    xlabel('hour')
end



% --- Executes on button press in pushbuttonLoadData.
function pushbuttonLoadData_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonLoadData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global figHandle2
LoadData
waitfor(figHandle2)


% --- Executes on selection change in popupmenuTempProfile.
function popupmenuTempProfile_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuTempProfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global update Plant Model_dir
update =0;
val = get(handles.popupmenuTempProfile,'value');
climate = get(handles.popupmenuTempProfile,'string');
filename = fullfile(Model_dir, 'component library','Weather',strcat(char(climate(val)),'.mat'));
load(filename)
climateName = {'1A'; '2A'; '2B'; '3A'; '3B'; '3B-Coast'; '3C'; '4A'; '4B'; '4C'; '5A'; '5B'; '6A'; '6B'; '7'; '8';};
if nnz(strcmp(climateName,char(climate(val))))>0
    Temperature = (Temperature/10-32)*5/9; %convert temperature data to celcius
end
if length(Temperature)~=length(Plant.Data.Timestamp)
    Steps = round(1/(Plant.Data.Timestamp(2)-Plant.Data.Timestamp(1)));%points in 1 day of data
    n = length(Plant.Data.Temperature);
    r = Steps/24; %data points per hour
    D = datevec(Plant.Data.Timestamp(1));
    Xs = 1+round(24*(Plant.Data.Timestamp(1)-datenum([D(1) 1 1])));%hour since jan1 of 1st point
    nS = floor(24*(Plant.Data.Timestamp(end)-Plant.Data.Timestamp(1)));%hours in dataset
    for t = 1:1:nS
        Plant.Data.Temperature(r*(t-1)+1:r*t) = Temperature(Xs);
        Xs = Xs+1;
        if Xs>length(Temperature)
            Xs=1;
        end
    end
    Plant.Data.Temperature(r*t:end) = Plant.Data.Temperature(r*t);
else
    Plant.Data.Temperature = Temperature;
end
Plant.Data.Weather=climate(val);
popupmenuAxes_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function popupmenuTempProfile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuTempProfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbuttonCancel.
function pushbuttonCancel_Callback(hObject, eventdata, ~)
% hObject    handle to pushbuttonCancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Plant figHandle DataOriginal
Plant.Data = DataOriginal;
close(figHandle)

% check if the day you are forecasting is during the week (mon through fri) and that it is not a holiday
function workday = isbusday(dateDay, Holidays, weekend)
weekendday = nnz(weekend);
    if weekday(dateDay)==weekendday
         dateDay=0;
    end
    for i=1:1:length(Holidays)
        if dateDay==Holidays(i)
            dateDay=0;
        end
    end
    if dateDay~=0
        dateDay=1;
    end
    workday=dateDay;
