function varargout = ViewResults2(varargin)
% VIEWRESULTS2 MATLAB code for ViewResults2.fig
%      VIEWRESULTS2, by itself, creates a new VIEWRESULTS2 or raises the existing
%      singleton*.
%
%      H = VIEWRESULTS2 returns the handle to a new VIEWRESULTS2 or the handle to
%      the existing singleton*.
%
%      VIEWRESULTS2('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VIEWRESULTS2.M with the given input arguments.
%
%      VIEWRESULTS2('Property','Value',...) creates a new VIEWRESULTS2 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ViewResults2_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ViewResults2_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ViewResults2RF

% Last Modified by GUIDE v2.5 07-Apr-2014 11:49:17

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ViewResults2_OpeningFcn, ...
                   'gui_OutputFcn',  @ViewResults2_OutputFcn, ...
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


% --- Executes just before ViewResults2 is made visible.
function ViewResults2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ViewResults2 (see VARARGIN)
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

global Model_dir Project
%plot options
plots1={'Electric Demand'
        'Electric Demand Hist'
        'Electric without HVAC'
        'Electric Demand by Day'
        'Cooling Demand'
        'Cooling Hist'
        'Heat Demand'
        'Heat Demand Hist'};
    
plots2 = {'NPC'; 'CO2'; 'NOx'; 'SO2';};
set(handles.popupmenuAxes1,'string',plots1,'value',1)
set(handles.popupmenuAxes2,'string',plots2,'value',1)

listUser = {'Residential'; 'Commercial';'Industrial';};
set(handles.popupmenuUser,'string',listUser,'value',2)
set(handles.checkboxScale,'value',1)
set(handles.checkboxGeneration,'value',1)

stateName = {'Alabama';'Alaska';'Arizona';'Arkansas';'California';'Colorado';'Connecticut';'Delaware';'Florida';'Georgia';
             'Hawaii';'Idaho';'Illinois';'Indiana';'Iowa';'Kansas';'Kentucky';'Louisiana';'Maine';'Maryland';
             'Massachusetts';'Michigan';'Minnesota';'Mississippi';'Missouri';'Montana';'Nebraska';'Nevada';'NewHampshire';'NewJersey';
             'NewMexico';'NewYork';'NorthCarolina';'NorthDakota';'Ohio';'Oklahoma';'Oregon';'Pennsylvania';'RhodeIsland';'SouthCarolina';
             'SouthDakota';'Tennessee';'Texas';'Utah';'Vermont';'Virginia';'Washington';'WestVirginia';'Wisconsin';'Wyoming';};
currentState = find(strcmp(Project.State,stateName));
set(handles.popupmenuState, 'string', stateName(:,1), 'value',currentState)

BuildType = {'Restaurant: full-service (sit down)';'Restaurant: quick-service (fast food)';'School: primary school';'School: secondary school';'Office: large office';'Office: medium office';'Office: small office';'Mid-rise apartment building';'Hospitality: large hotel';'Hospitality: small hotel/motel';'Health care: large hospital';'Health care: outpatient facility';'Retail: big-box, standalone retail store';'Retail: retail store located in a strip mall';'Retail: supermarket';'Unrefrigerated warehouse';'Multiple';'Real Building'};
BuildList = {'SDRest'; 'FFRest'; 'Sch-pri'; 'Sch-sec'; 'LgOff'; 'MdOff'; 'SmOff'; 'MRapt'; 'LgHotel'; 'SmHotel'; 'Hospital'; 'OutP'; 'Retail'; 'StMall'; 'SMarket'; 'ware';};
if length(Project.Building.NameList)>1
    currentBuild = 17;
else
    name = char(Project.Building.NameList{1});
    ind = strfind(name,'_');
    currentBuild = find(strcmp(name(1:ind(1)-1),BuildList),1);
    if isempty(currentBuild)
        currentBuild = 18;
    end
end
set(handles.popupmenuBuilding,'string',BuildType,'value',currentBuild)

Vintage = {'2010 construction (ASHRAE 90.1-2010)';'2007 construction (ASHRAE 90.1-2007)';'2004 construction 90.1-2004';'“Post-1980” construction (ASHRAE 90.1-1989)';'“Pre-1980” construction';'Real Building';};
VintageList = {'New2010';'New2007';'New2004';'Post1980';'Pre1980';};
if currentBuild == 18;
    currentVintage = 6;
else
    name = char(Project.Building.NameList{1});
    ind = strfind(name,'_');
    currentVintage = find(strcmp(name(ind(2)+1:end),VintageList),1);
    if isempty(currentVintage)
        currentVintage = 1;
    end
end
set(handles.popupmenuVintage,'string',Vintage,'value',currentVintage)

Control_dir=fullfile(Model_dir, 'System Library', 'Control'); 
ControlFiles=dir(fullfile(Control_dir,'*.mat'));
ControlList=strrep({ControlFiles.name},'.mat','');
currentControl = find(strcmp(Project.Control.Name,ControlList));
if isempty(currentControl)
    currentControl = 1;
end
set(handles.popupmenuControl,'string',ControlList,'value',currentControl)

Grid_dir=fullfile(Model_dir,'System Library','Grid'); 
GridFiles=dir(fullfile(Grid_dir,'*.mat'));
GridList=strrep({GridFiles.name},'.mat','');
currentUtility = find(strcmp(Project.Utilities.Name,GridList));
if isempty(currentUtility)
    currentUtility = 1;
end
set(handles.popupmenuUtility,'string',GridList,'value',currentUtility)

set(handles.sliderZoom,'Min',1,'Max',4,'Value',1,'SliderStep',[1/3,1/3])
sliderZoom_Callback(handles.sliderZoom, eventdata, handles)

BoilerFuel = {'NaturalGas','Oil','Coal','Propane'};
set(handles.popupmenuBoilerType,'string',BoilerFuel,'value',1)

set(handles.editSellBack,'string',num2str(Project.Utilities.Grid.SellBackRate))
if ~isfield(Project.Utilities.Grid,'SellBackRatePerc');
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

updateViewResults(handles,0)


% --- Outputs from this function are returned to the command line.
function varargout = ViewResults2_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in popupmenuAxes1.
function popupmenuAxes1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuAxes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

val=get(hObject,'value');
str=get(hObject,'string');
currentstr=str{val};

global Project
Ts = 8760/length(Project.Building.DemandE);
axes(handles.axes1)
cla reset
set(gca,'tag','axes1') 

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
    if StartDay+6<=monthDays(StartMonth+1)
        plotAxis = datenum([2014*ones(8,1) StartMonth*ones(8,1)  (StartDay:StartDay+7)'  0*ones(8,1)  0*ones(8,1) 0*ones(8,1)]);
        xlabel(monthLabel(StartMonth,:))
    else next = StartDay+6-monthDays(StartMonth+1);
        plotAxis = datenum([2014*ones(8,1) StartMonth*ones(8,1)  [(StartDay:month(StartMonth+1)) 1:next+1]'  0*ones(8,1)  0*ones(8,1) 0*ones(8,1)]);
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
        PlotType = 'plot';
        Y = Project.Building.DemandH(1+24/Ts*(day1-1):24/Ts*lastDay);
    case 'Heat Demand Hist'
        PlotType = 'hist';
        Y = Project.Building.DemandH(1+24/Ts*(day1-1):24/Ts*lastDay);
    case 'Electric Demand'
        PlotType = 'plot';
        Y = Project.Building.DemandE(1+24/Ts*(day1-1):24/Ts*lastDay);
    case 'Electric Demand by Day'
        PlotType = 'pcolor';
        xgrid=floor(Project.Building.Hour(1+24/Ts*(day1-1))):floor(Project.Building.Hour(24/Ts*lastDay));
        ygrid=Ts:Ts:24;
        z=NaN*ones(length(xgrid),length(ygrid));
        for i=1:lastDay-day1+1
            z(i,:)=Project.Building.DemandE(1+24/Ts*(i-1+day1-1):(i+day1-1)*24/Ts);
        end
    case 'Electric Demand Hist'
        PlotType = 'hist';
        Y = Project.Building.DemandE(1+24/Ts*(day1-1):24/Ts*lastDay);
    case 'Cooling Demand'
        PlotType = 'plot';
        Y = Project.Building.DemandC(1+24/Ts*(day1-1):24/Ts*lastDay);
    case 'Cooling Hist'
        PlotType = 'hist';
        Y = Project.Building.DemandC(1+24/Ts*(day1-1):24/Ts*lastDay);
    case 'Electric without HVAC'
        PlotType = 'plot';
        Y = Project.Building.DemandE(1+24/Ts*(day1-1):24/Ts*lastDay)-Project.Building.CoolingElectricalLoad(1+24/Ts*(day1-1):24/Ts*lastDay);
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
    h = colorbar;
    xlabel(h,'(kW)');
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
elec1 = strcmp(str(get(handles.popupmenuAxes1,'value')),'Electric Demand') || strcmp(str(get(handles.popupmenuAxes1,'value')),'Electric without HVAC');
if elec1
    set(handles.checkboxGeneration,'enable','on')
else set(handles.checkboxGeneration,'enable','off')
end
if isfield(Project.System,'TES')
    set(handles.checkboxTESshift,'enable','on')
else set(handles.checkboxTESshift,'enable','off')
end
if strcmp(str(get(hObject,'value')),'Electric Demand') || strcmp(str(get(hObject,'value')),'Electric without HVAC')
    if get(handles.checkboxGeneration,'Value')==1
        genPow = Project.Result.eOut.GenPower(1+24/Ts*(day1-1):24/Ts*lastDay);
        hold on
        plot(X,genPow,'r');
        hold off
    end
    if get(handles.checkboxTESshift,'Value')==1 && isfield(Project.System,'TES')
        hold on
        plot(X,Project.Result.eOut.DemandEshift(1+24/Ts*(day1-1):24/Ts*lastDay),'g');
        hold off
    end
end

% --- Executes during object creation, after setting all properties.
function popupmenuAxes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuAxes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
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
popupmenuAxes1_Callback(handles.popupmenuAxes1,eventdata,handles)

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
popupmenuAxes1_Callback(handles.popupmenuAxes1,eventdata,handles)

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
popupmenuAxes1_Callback(handles.popupmenuAxes1,eventdata,handles)

% --- Executes on button press in checkboxTESshift.
function checkboxTESshift_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxTESshift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
popupmenuAxes1_Callback(handles.popupmenuAxes1,eventdata,handles)


% --- Executes on selection change in popupmenuAxes1.
function popupmenuAxes2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuAxes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
val=get(handles.popupmenuAxes2,'value');
str=get(handles.popupmenuAxes2,'string');
currentstr=str{val};
global RESULTS 
axes(handles.axes2);

StateNum = get(handles.popupmenuState,'Value');
stateAbrev = get(handles.popupmenuState,'string');
State = stateAbrev(StateNum);
A = get(handles.uipanelGridEmissions,'SelectedObject');
switch get(A,'Tag')
    case 'radiobuttonStateAvg'
        GridMix = 1;
    case 'radiobuttonCombustionOnly'
        GridMix = 2;
end
Boiler = get(handles.popupmenuBoilerType,'value');
[BaselineEmission DispatchEmission] = EmissionsCalculated(RESULTS,State,GridMix,Boiler);

       
switch currentstr
case 'NPC'
        bar([RESULTS.costOut.NPVbaselineDemCharges RESULTS.costOut.NPVbaselineUseCharges  RESULTS.costOut.NPVbaselineFuelCost    RESULTS.costOut.NPVbaselineOandM RESULTS.costOut.NPVbaselineFinance;...
             RESULTS.costOut.NPVnewDemCharges RESULTS.costOut.NPVnewUseCharges RESULTS.costOut.NPVnewFuelCost  RESULTS.costOut.NPVnewOandM RESULTS.costOut.NPVnewFinance],'stacked')
        ylabel('Fuel, Grid, and O&M/Finance Costs ($)')
        set(gca,'XTickLabel',{'Baseline','System'})
        xlim([0.5 2.5])
        legend('Demand', 'Grid' ,'Fuel', 'O&M','Financing','Location','NorthEastOutside')
	case 'CO2'
        bar([BaselineEmission(1,1:3);DispatchEmission(1,1:3);],'stacked')
        ylabel('Emissions of CO2 (tons)')
        set(gca,'XTickLabel',{'Baseline','System'})
        xlim([0.5 2.5])
        legend('Grid', 'CHP' ,'Boiler', 'Location','NorthEastOutside')
    case 'NOx'
        bar([BaselineEmission(2,1:3);DispatchEmission(2,1:3);],'stacked')
        ylabel('Emissions of  NOx (lbs)')
        set(gca,'XTickLabel',{'Baseline','System'})
        xlim([0.5 2.5])
        legend('Grid', 'CHP' ,'Boiler', 'Location','NorthEastOutside')
    case 'SO2'
        bar([BaselineEmission(3,1:3);DispatchEmission(3,1:3);],'stacked')
        ylabel('Emissions of SO2 (lbs)')
        set(gca,'XTickLabel',{'Baseline','System'})
        xlim([0.5 2.5])
        legend('Grid', 'CHP' ,'Boiler', 'Location','NorthEastOutside')
end

% --- Executes during object creation, after setting all properties.
function popupmenuAxes2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuAxes2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function sliderCost_Callback(hObject, eventdata, handles)
% hObject    handle to sliderCost (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Project
Project.Economic.InstallCost=round(get(hObject,'Value'));
updateViewResults(handles,1)

% --- Executes during object creation, after setting all properties.
function sliderCost_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderCost (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function editCost_Callback(hObject, eventdata, handles)
% hObject    handle to editCost (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Project
installCost=str2double(get(hObject,'string'));
Project.Economic.InstallCost=installCost;
updateViewResults(handles,1)

% --- Executes during object creation, after setting all properties.
function editCost_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editCost (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function sliderIncent_Callback(hObject, eventdata, handles)
% hObject    handle to sliderIncent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Project
Project.Economic.Incentive=round(get(hObject,'Value'));
updateViewResults(handles,1)

% --- Executes during object creation, after setting all properties.
function sliderIncent_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderIncent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function editIncent_Callback(hObject, eventdata, handles)
% hObject    handle to editIncent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Project
Project.Economic.Incentive=str2double(get(hObject,'string'));
updateViewResults(handles,1)

% --- Executes during object creation, after setting all properties.
function editIncent_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editIncent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function sliderLifespan_Callback(hObject, eventdata, handles)
% hObject    handle to sliderLifespan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Project
Project.Economic.LifeYrs=round(get(hObject,'Value'));
CHPSize = 0;
for i = 1:1:length(Project.System.CHP)
    CHPSize = CHPSize+Project.System.CHP(i).SysSize(1);
end
Project.Economic.LifekWh = Project.Economic.LifeYrs*CHPSize*8760;
updateViewResults(handles,1)

% --- Executes during object creation, after setting all properties.
function sliderLifespan_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderLifespan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


function editLifespan_Callback(hObject, eventdata, handles)
% hObject    handle to editLifespan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Project
Project.Economic.LifeYrs=str2double(get(hObject,'string'));
CHPSize = 0;
for i = 1:1:length(Project.System.CHP)
    CHPSize = CHPSize+Project.System.CHP(i).SysSize(1);
end
Project.Economic.LifekWh = Project.Economic.LifeYrs*CHPSize*8760;
updateViewResults(handles,1)

% --- Executes during object creation, after setting all properties.
function editLifespan_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editLifespan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function sliderPayback_Callback(hObject, eventdata, handles)
% hObject    handle to sliderPayback (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Project
Project.Economic.Payback=round(get(hObject,'Value'));
updateViewResults(handles,1)

% --- Executes during object creation, after setting all properties.
function sliderPayback_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderPayback (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function editPayback_Callback(hObject, eventdata, handles)
% hObject    handle to editPayback (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Project
Project.Economic.Payback=str2double(get(hObject,'string'));
OriginCost = Project.Economic.InstallCost;
tempPayback = calcPayback();
error = Project.Economic.Payback-tempPayback;
while abs(error)>.1
    if tempPayback >=19 || Project.Economic.InstallCost<50
        error = 0;
        Project.Economic.Payback = 20;
        Project.Economic.InstallCost = OriginCost;
    else
        Project.Economic.InstallCost = Project.Economic.InstallCost*((1+error/tempPayback)^.5);
        tempPayback = calcPayback();
        error = Project.Economic.Payback-tempPayback;
    end
end
updateViewResults(handles,1)

% --- Executes during object creation, after setting all properties.
function editPayback_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editPayback (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on slider movement.
function sliderCHPsize_Callback(hObject, eventdata, handles)
% hObject    handle to sliderCHPsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Project
CHPSize = 0;
for i = 1:1:length(Project.System.CHP)
    CHPSize = CHPSize+Project.System.CHP(i).SysSize(1);
end
NewSize = round(get(hObject,'Value'));
Project.Economic.LifekWh = Project.Economic.LifeYrs*NewSize*8760;
for i = 1:1:length(Project.System.CHP)
    Project.System.CHP(i).SysSize(1)= Project.System.CHP(i).SysSize(1)*NewSize/CHPSize;
end
updateViewResults(handles,0)

% --- Executes during object creation, after setting all properties.
function sliderCHPsize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderCHPsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function editCHPsize_Callback(hObject, eventdata, handles)
% hObject    handle to editCHPsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Project
NewSize = str2double(get(hObject,'string'));
CHPSize = 0;
for i = 1:1:length(Project.System.CHP)
    CHPSize = CHPSize+Project.System.CHP(i).SysSize(1);
end
Project.Economic.LifekWh = Project.Economic.LifeYrs*NewSize*8760;
for i = 1:1:length(Project.System.CHP)
    Project.System.CHP(i).SysSize(1)= Project.System.CHP(i).SysSize(1)*NewSize/CHPSize;
    Project.System.CHP(i).SysSize(2) = Project.System.CHP(i).SysSize(1)/Project.System.CHP(i).TurnDown;
end
updateViewResults(handles,0)

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
MSG = msgbox('Re-sizing CHP System');
ScaleSystemComponents(100,100,100,120,1,2,1,2,2)
updateViewResults(handles,0)
close(MSG)

% --- Executes on selection change in popupmenuBuilding.
function popupmenuBuilding_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuBuilding (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Project Model_dir
type=get(handles.popupmenuBuilding, 'value');
vintage=get(handles.popupmenuVintage, 'value');

StateNum = get(handles.popupmenuState,'Value');
stateAbrev = {'AL';'AK';'AZ';'AR';'CA';'CO';'CT';'DE';'FL';'GA';'HI';'ID';'IL';'IN';'IA';'KS';'KY';'LA';'ME';'MD';'MA';'MI';'MN';'MS';'MO';
              'MT';'NE';'NV';'NH';'NJ';'NM';'NY';'NC';'ND';'OH';'OK';'OR';'PA';'RI';'SC';'SD';'TN';'TX';'UT';'VT';'VA';'WA';'WV';'WI';'WY';};
State = stateAbrev(StateNum);
ASHRAEZone = ClimateZone( char(State),00000);
Climates = {'1A'; '2A'; '2B'; '3A'; '3B'; '3B-Coast'; '3C'; '4A'; '4B'; '4C'; '5A'; '5B'; '6A'; '6B'; '7'; '8';};
climate = find(strcmp(ASHRAEZone,Climates),1);

vintageName = {'New2010';'New2007';'New2004';'Post1980';'Pre1980';};
climateName = {'_1A_'; '_2A_'; '_2B_'; '_3A_'; '_3B_'; '_3B-Coast_'; '_3C_'; '_4A_'; '_4B_'; '_4C_'; '_5A_'; '_5B_'; '_6A_'; '_6B_'; '_7_'; '_8_';};
buildTypeName = {'SDRest'; 'FFRest'; 'Sch-pri'; 'Sch-sec'; 'LgOff'; 'MdOff'; 'SmOff'; 'MRapt'; 'LgHotel'; 'SmHotel'; 'Hospital'; 'OutP'; 'Retail'; 'StMall'; 'SMarket'; 'ware';};
buildDir = fullfile(Model_dir,'System Library', 'Buildings');
buildCHPtemp = Project.Building.CHPtemp;
Project.Building.DemandE = zeros(35040,1);
Project.Building.DemandC = zeros(35040,1);
Project.Building.DemandH = zeros(35040,1);
Project.Building.CoolingElectricalLoad = zeros(35040,1);
for i = 1:1:length(Project.Building.NameList)
    name = char(Project.Building.NameList{i});
    ind =strfind(name,'_');
    oldClimate = name(ind(1):ind(2));
    oldVintage = name(ind(2)+1:end);
    oldBuild = name(1:ind(1)-1);
    name = char(strrep(name,oldBuild,buildTypeName(type)));
    name = char(strrep(name,oldClimate,climateName(climate)));
    name = char(strrep(name,oldVintage,vintageName(vintage)));
    load(strcat(buildDir,filesep,name))
    component.DistHeat = get(handles.checkboxDistHeat,'value');
    component.DistCool = get(handles.checkboxDistCool,'value');
    Project.Building.DemandE = Project.Building.DemandE+ component.DemandE;
    Project.Building.DemandC = Project.Building.DemandC+ component.DistCool*component.DemandC;
    Project.Building.DemandH = Project.Building.DemandH+ component.DistHeat*component.DemandH;
    Project.Building.NonDistH = Project.Building.NonDistH + (1-component.DistHeat)*component.DemandH;
    Project.Building.CoolingElectricalLoad = Project.Building.CoolingElectricalLoad+ component.DistCool*component.CoolingElectricalLoad;
    Project.Building.NameList{i} = cellstr(name);
    Project.Building.CHPtemp(i) = buildCHPtemp(i);
end
updateViewResults(handles,0)

% --- Executes during object creation, after setting all properties.
function popupmenuBuilding_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuBuilding (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in popupmenuVintage.
function popupmenuVintage_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuVintage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
popupmenuBuilding_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function popupmenuVintage_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuVintage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in checkboxDistCool.
function checkboxDistCool_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxDistCool (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
popupmenuBuilding_Callback(hObject, eventdata, handles)

% --- Executes on button press in checkboxDistHeat.
function checkboxDistHeat_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxDistHeat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
popupmenuBuilding_Callback(hObject, eventdata, handles)


% --- Executes on selection change in popupmenuControl.
function popupmenuControl_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuControl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Project Model_dir
val=get(hObject,'value');
str=get(hObject,'string');
currentstr=str{val};
Control_dir=fullfile(Model_dir,'System Library','Control'); %strrep(which('NREL_FCModel.m'),'\main\NREL_FCModel.m','\component library\Control');
load(strcat(Control_dir,filesep,currentstr))
Project.Control = component;
updateViewResults(handles,0)

% --- Executes during object creation, after setting all properties.
function popupmenuControl_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuControl (see GCBO)
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
global Project  Model_dir
val=get(hObject,'value');
str=get(hObject,'string');
currentstr=str{val};
Grid_dir=fullfile(Model_dir,'System Library','Grid'); % strrep(which('NREL_FCModel.m'),'\main\NREL_FCModel.m','\component library\Grid');
load(strcat(Grid_dir,filesep,currentstr))
Project.Utilities.Grid = component;
updateViewResults(handles,0)

% --- Executes during object creation, after setting all properties.
function popupmenuUtility_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuUtility (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in popupmenuBoilerType.
function popupmenuBoilerType_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuBoilerType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
updateViewResults(handles,1)

% --- Executes during object creation, after setting all properties.
function popupmenuBoilerType_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuBoilerType (see GCBO)
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
global Project
list = get(handles.popupmenuState,'string');
val = get(handles.popupmenuState,'value');
Project.State = list(val);
popupmenuBuilding_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function popupmenuState_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuState (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in checkboxScale.
function checkboxScale_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxScale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
updateViewResults(handles,1)

% --- Executes on selection change in popupmenuUser.
function popupmenuUser_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuUser (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
updateViewResults(handles,1)

% --- Executes during object creation, after setting all properties.
function popupmenuUser_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuUser (see GCBO)
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
set(handles.uipanelSellBack,'SelectedObject',handles.radiobuttonFixedSellBack)
updateViewResults(handles,1)

% --- Executes during object creation, after setting all properties.
function editSellBack_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editSellBack (see GCBO)
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
updateViewResults(handles,1)


% --- Executes when selected object is changed in uipanelChillType.
function uipanelGridEmissions_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanelChillType 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
updateViewResults(handles,1)

%% 
function updateViewResults(handles,skipRun)
global Project RESULTS
MSGvr = msgbox('Recalculating');

str = get(handles.popupmenuUser,'string');
val = get(handles.popupmenuUser,'value');
user = str{val};
scale = get(handles.checkboxScale,'value');
if skipRun
    Project.Result.costOut=FinancialCalcs2(Project.Result.Baseline,Project.Result.Dispatch,Project.Utilities.Grid,Project.Economic,Project.State,user,Project.Utilities.NatGas,scale);
else
    Project.Result = runAnalyses(Project,user,scale);
end

RESULTS = Project.Result;
eventdata = 0;
popupmenuAxes1_Callback(handles.popupmenuAxes1,eventdata,handles)
popupmenuAxes2_Callback(handles.popupmenuAxes2,eventdata,handles)

BaseGridCost = RESULTS.costOut.Baseline.TotalGridBill/sum(RESULTS.Baseline.Elec);
CaseGridCost = RESULTS.costOut.Dispatch.TotalGridBill/sum(RESULTS.Dispatch.Elec);
SelfGenValue = (RESULTS.costOut.Baseline.TotalGridBill-RESULTS.costOut.Dispatch.TotalGridBill)/RESULTS.Dispatch.ElecTotProd;

CHPSize = 0;
for i = 1:1:length(Project.System.CHP)
    CHPSize = CHPSize+Project.System.CHP(i).SysSize(1);
end

set(handles.editCost,'string',round(Project.Economic.InstallCost))
set(handles.editIncent,'string',round(Project.Economic.Incentive))
set(handles.editLifespan,'string',round(Project.Economic.LifeYrs))
set(handles.editPayback,'string',round(RESULTS.costOut.Payback))
set(handles.editCHPsize,'string',round(CHPSize))
set(handles.sliderCost,'Min',0,'Max',12000,'Value',Project.Economic.InstallCost,'SliderStep',[1/24,1/24])
set(handles.sliderIncent,'Min',0,'Max',Project.Economic.InstallCost,'Value',Project.Economic.Incentive,'SliderStep',[1/20,1/20])
set(handles.sliderLifespan,'Min',1,'Max',20,'Value',Project.Economic.LifeYrs,'SliderStep',[1/20,1/20])
set(handles.sliderPayback,'Min',0,'Max',20,'Value',RESULTS.costOut.Payback,'SliderStep',[1/20,1/20])
set(handles.sliderCHPsize,'Min',0,'Max',2*CHPSize+1,'Value',CHPSize,'SliderStep',[1/20,1/20])
Project.Economic.Payback = RESULTS.costOut.Payback;

cost = Project.Result.costOut;
costNum = [mean(cost.BaselineGridBill+cost.BaselineFuel);mean(cost.CostPerYear);mean(cost.FuelCosts);mean(cost.NewGridBill);];
if min(costNum)<1000
    costNum = round(costNum);
elseif min(costNum)<1e4
    costNum = 10*round(costNum/10);
elseif min(costNum)<1e5
    costNum = 100*round(costNum/100);
elseif min(costNum)<1e6
    costNum = 1e3*round(costNum/1e3);
elseif min(costNum)<1e7
    costNum = 1e4*round(costNum/1e4);
else costNum = 1e5*round(costNum/1e5);
end

ResultsLabel = {strcat('CHP Size (kW):');
    strcat('CHP Capacity Factor:');
    strcat('Self GenerationPproportion: ');
    strcat('Proportion of Heating Met by CHP:');
    strcat('Annual cost, Baseline / with CHP ($):');
    strcat('CHP Annual Cost Electric Utility/CHP fuel ($)');
    strcat('Avg Grid Cost, Baseline/ with CHP (cents/kWh):');
    strcat('Avg CHP Electricity Cost (cents/kWh):');};
Results = {strcat(num2str(round(RESULTS.eOut.SysSize)));
    strcat(num2str(RESULTS.eOut.Total_Electricity_produced_kWh/(RESULTS.eOut.SysSize*8760),2));
    strcat(num2str(RESULTS.eOut.SelfGen,2));
    strcat(num2str(RESULTS.eOut.Total_DG_Heat_out_kWh/(RESULTS.eOut.Total_DG_Heat_out_kWh+ RESULTS.eOut.Total_peak_burner_heat_out_kWh),2));
    strcat(num2str(round(costNum(1))),' / ',num2str(round(costNum(2))));
    strcat(num2str(round(costNum(4))),' / ',num2str(round(costNum(3))));
    strcat(num2str(BaseGridCost*100,3),' / ',num2str(CaseGridCost*100,3));
    strcat(num2str(SelfGenValue*100,3));};

set(handles.textResults,'string',Results)
set(handles.textResultsLabel,'string',ResultsLabel)

vintageName = {'New2010';'New2007';'New2004';'Post1980';'Pre1980';};
climateName = {'_1A_'; '_2A_'; '_2B_'; '_3A_'; '_3B_'; '_3B-Coast_'; '_3C_'; '_4A_'; '_4B_'; '_4C_'; '_5A_'; '_5B_'; '_6A_'; '_6B_'; '_7_'; '_8_';};
buildTypeName = {'SDRest'; 'FFRest'; 'Sch-pri'; 'Sch-sec'; 'LgOff'; 'MdOff'; 'SmOff'; 'MRapt'; 'LgHotel'; 'SmHotel'; 'Hospital'; 'OutP'; 'Retail'; 'StMall'; 'SMarket'; 'ware';};

if length(Project.Building.NameList)>1
    BuildVal(1) = 17;
    name = char(Project.Building.NameList{1});
    ind =strfind(name,'_');
    BuildVal(2) = find(strcmp(name(ind(1):ind(2)),climateName));
    BuildVal(3) = find(strcmp(name(ind(2)+1:end),vintageName));
else name = char(Project.Building.NameList{1});
    ind =strfind(name,'_');
    if isempty(find(strcmp(name(1:ind(1)-1),buildTypeName),1));
        BuildVal(1:3) = [18, 16, 3];
    else
        BuildVal(1) = find(strcmp(name(1:ind(1)-1),buildTypeName));
        BuildVal(2) = find(strcmp(name(ind(1):ind(2)),climateName));
        BuildVal(3) = find(strcmp(name(ind(2)+1:end),vintageName));
    end
end

set(handles.popupmenuBuilding,'value',BuildVal(1))
set(handles.popupmenuVintage,'value',BuildVal(3))
set(handles.checkboxDistHeat,'value',Project.Building.DistHeat(1))
set(handles.checkboxDistCool,'value',Project.Building.DistCool(1))

ControlList = get(handles.popupmenuControl,'string');
GridList = get(handles.popupmenuUtility,'string');
ControlVal = find(strcmp(Project.Control.Name,ControlList));
GridVal = find(strcmp(Project.Utilities.Grid.Name,GridList));

set(handles.popupmenuControl,'value',ControlVal)
set(handles.popupmenuUtility,'value',GridVal)
close(MSGvr)


function payBack = calcPayback()
global Project
Project.Result = runAnalyses(Project);
payBack = Project.Result.costOut.Payback;


% --- Executes on button press in pushbuttonEmissionDetail.
function pushbuttonEmissionDetail_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonEmissionDetail (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
EmissionControl()


function editSellbackPerc_Callback(hObject, eventdata, handles)
% hObject    handle to editSellbackPerc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Project
Project.Utilities.Grid.SellBackPerc =str2double(get(hObject,'string'));
updateViewResults(handles,1)

% --- Executes during object creation, after setting all properties.
function editSellbackPerc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editSellbackPerc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbuttonExport.
function pushbuttonExport_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonExport (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Boiler = get(handles.popupmenuBoilerType,'value');
ExportResult('ViewResults',Boiler)


% --- Executes on button press in pushbuttonGHGsize.
function pushbuttonGHGsize_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonGHGsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
list={'CO2'; 'SO2'; 'NOx';};
    [selection, ok]=listdlg('ListString',list);
    if ok
        selection=list{selection(1)};
        switch selection
            case 'CO2'
                autoSizing('CO2','CHP')
            case 'SO2'
                autoSizing('SO2','CHP')
            case 'NOx'
                autoSizing('NOx','CHP')
        end
    end
updateViewResults(handles,0)
    
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
updateViewResults(handles,0)
