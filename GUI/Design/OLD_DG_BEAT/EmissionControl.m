function varargout = EmissionControl(varargin)
% EMISSIONCONTROL MATLAB code for EmissionControl.fig
%      EMISSIONCONTROL, by itself, creates a new EMISSIONCONTROL or raises the existing
%      singleton*.
%
%      H = EMISSIONCONTROL returns the handle to a new EMISSIONCONTROL or the handle to
%      the existing singleton*.
%
%      EMISSIONCONTROL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EMISSIONCONTROL.M with the given input arguments.
%
%      EMISSIONCONTROL('Property','Value',...) creates a new EMISSIONCONTROL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before EmissionControl_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to EmissionControl_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ViewResults2RF

% Last Modified by GUIDE v2.5 11-Feb-2014 16:04:34

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @EmissionControl_OpeningFcn, ...
                   'gui_OutputFcn',  @EmissionControl_OutputFcn, ...
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


% --- Executes just before EmissionControl is made visible.
function EmissionControl_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to EmissionControl (see VARARGIN)

% Choose default command line output for EmissionControl
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

global Model_dir Project

plots1={'Electric'
        'CO2 Emissions'
        'SO2 Emissions'
        'NOx Emissions'
        'Grid emission factors'};
    
set(handles.popupmenuAxes1,'string',plots1,'value',2)

stateName = {'Alabama';'Alaska';'Arizona';'Arkansas';'California';'Colorado';'Connecticut';'Delaware';'Florida';'Georgia';
             'Hawaii';'Idaho';'Illinois';'Indiana';'Iowa';'Kansas';'Kentucky';'Louisiana';'Maine';'Maryland';
             'Massachusetts';'Michigan';'Minnesota';'Mississippi';'Missouri';'Montana';'Nebraska';'Nevada';'NewHampshire';'NewJersey';
             'NewMexico';'NewYork';'NorthCarolina';'NorthDakota';'Ohio';'Oklahoma';'Oregon';'Pennsylvania';'RhodeIsland';'SouthCarolina';
             'SouthDakota';'Tennessee';'Texas';'Utah';'Vermont';'Virginia';'Washington';'WestVirginia';'Wisconsin';'Wyoming';};
set(handles.popupmenuState, 'string', stateName(:,1), 'value',5)


Control_dir=fullfile(Model_dir, 'System Library', 'Control'); 
ControlFiles=dir(fullfile(Control_dir,'*.mat'));
ControlList=strrep({ControlFiles.name},'.mat','');
currentControl = find(strcmp(Project.Control.Name,ControlList));
if isempty(currentControl)
    currentControl = 1;
end
set(handles.popupmenuControl,'string',ControlList,'value',currentControl)

set(handles.sliderZoom,'Min',1,'Max',4,'Value',1)
set(handles.sliderDate,'Min',1,'Max',1.1,'Value',1)
% sliderZoom_Callback(handles.sliderZoom, eventdata, handles)

plotSubregion={'AKGD';'AKMS';'ERCT';'FRCC';'HIMS';'HIOA';'MROE';'MROW';'NYLI';'NEWE';'NYCW';'NYUP';'RFCE';'RFCM';'RFCW';'SRMW';'SRMV';'SRSO';'SRTV';'SRVC';'SPNO';'SPSO';'CAMX';'NWPP';'RMPA';'AZNM'};
set(handles.popupmenuSubregion,'string',plotSubregion,'value',23)
BoilerFuel = {'NaturalGas','Oil','Coal','Propane'};
set(handles.popupmenuBoilerType,'string',BoilerFuel,'value',1)
updateViewResults(handles,1)


% --- Outputs from this function are returned to the command line.
function varargout = EmissionControl_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


function updateViewResults(handles,skipRun)
global Project RESULTS
MSG = msgbox('Recalculating');
if ~skipRun
    Project.Result = runAnalyses(Project);
end
RESULTS = Project.Result;
eventdata = 0;

ControlList = get(handles.popupmenuControl,'string');
ControlVal = find(strcmp(Project.Control.Name,ControlList));
set(handles.popupmenuControl,'value',ControlVal)

popupmenuAxes1_Callback(handles.popupmenuAxes1,eventdata,handles)

axes(handles.axes2);

Boiler = get(handles.popupmenuBoilerType,'value');

StateNum = get(handles.popupmenuState,'Value');
stateAbrev = {'AL';'AK';'AZ';'AR';'CA';'CO';'CT';'DE';'FL';'GA';'HI';'ID';'IL';'IN';'IA';'KS';'KY';'LA';'ME';'MD';'MA';'MI';'MN';'MS';'MO';
              'MT';'NE';'NV';'NH';'NJ';'NM';'NY';'NC';'ND';'OH';'OK';'OR';'PA';'RI';'SC';'SD';'TN';'TX';'UT';'VT';'VA';'WA';'WV';'WI';'WY';};
State = stateAbrev(StateNum);
[BaselineEmission DispatchEmission] = EmissionsCalculated(RESULTS,State,1,Boiler);

                
LifekWh = Project.Economic.LifekWh;
if isempty(LifekWh)
    CHPsize = 0;
    for n = 1:1:length(Project.System.CHP)
        CHPsize = CHPsize+Project.System.CHP(n).SysSize(1,1);
    end
    LifekWh = CHPsize*8760*Project.Economic.LifeYrs;
end
ActualYears = LifekWh/RESULTS.eOut.Total_Electricity_produced_kWh;
GroupedData = ActualYears*[BaselineEmission DispatchEmission];
bar(GroupedData,'grouped')
ylabel('Emissions of CO2,NOx,andSO2 (tons,lbs,lbs)')
set(gca,'XTickLabel',{'CO2','NOx','SO2'})
% xlim([0.5 2.5])
legend('Baseline', 'With CHP', 'Location','NorthEastOutside')

close(MSG)


% --- Executes on selection change in popupmenuAxes1.
function popupmenuAxes1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuAxes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

val=get(hObject,'value');
str=get(hObject,'string');
currentstr=str{val};

global Project RESULTS
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

Boiler = get(handles.popupmenuBoilerType,'value');
BoilerCO2 = [0.3988 0.5565 0.7066 0.4742];%lb CO2 per kWh
BoilerCO2 = BoilerCO2(Boiler);
BoilerNOx = 7.7e-4; %lbNOx/kWh
BoilerSO2 = 2.5e-3; %lb SO2/kWh

%retrive emission profile
DemandE = Project.Building.DemandE;
DemandH = Project.Building.DemandH;
BaselinePurchaseHour = zeros(1,8760);
BaselineHeatHour = zeros(1,8760);
for i = 1:1:8760
    BaselinePurchaseHour(i) = sum(DemandE(1+1/Ts*(i-1):1/Ts*i)*Ts);
    BaselineHeatHour(i) = sum(DemandH(1+1/Ts*(i-1):1/Ts*i)*Ts);
end
% eGridZones = get(handles.popupmenuSubregion,'string');
% val = get(handles.popupmenuSubregion,'value');
% eGridRegion = eGridZones(val);
StateNum = get(handles.popupmenuState,'Value');
stateAbrev = {'AL';'AK';'AZ';'AR';'CA';'CO';'CT';'DE';'FL';'GA';'HI';'ID';'IL';'IN';'IA';'KS';'KY';'LA';'ME';'MD';'MA';'MI';'MN';'MS';'MO';
              'MT';'NE';'NV';'NH';'NJ';'NM';'NY';'NC';'ND';'OH';'OK';'OR';'PA';'RI';'SC';'SD';'TN';'TX';'UT';'VT';'VA';'WA';'WV';'WI';'WY';};
State = stateAbrev(StateNum);
[CO2, NOx, SO2] = EmissionProfile(State);
clear X
X = datenum(2014,1,1,0,0,0)+((day1-1)+(Ts/24):(1/24):lastDay)';
switch currentstr
    case 'Electric'
        X = datenum(2014,1,1,0,0,0)+((day1-1)+(Ts/24):(Ts/24):lastDay)';
        Y=Project.Building.DemandE(1+24/Ts*(day1-1):24/Ts*lastDay);
        Y2= Project.Result.eOut.GenPower(1+24/Ts*(day1-1):24/Ts*lastDay);
        Ylab = 'Demand (kW)'; 
        LegendText = {'Demand';'CHP Generation';};
    case 'CO2 Emissions'
        baselineProfile = BaselinePurchaseHour.*CO2' + BoilerCO2*BaselineHeatHour;
        Y = baselineProfile(1+24*(day1-1):24*lastDay);
        annualProfile = Project.Result.eOut.Grid_purchases_hourly_kWh.*CO2 + RESULTS.eOut.CO2 + BoilerCO2*RESULTS.eOut.BoilerHeatHour;
        Y2 = annualProfile(1+24*(day1-1):24*lastDay);
        Ylab = 'CO2 Emissions (tons/hr)';
        LegendText = {'Baseline';'With CHP';};
    case 'SO2 Emissions'
        baselineProfile = BaselinePurchaseHour.*SO2' + BoilerSO2*BaselineHeatHour;
        Y = baselineProfile(1+24*(day1-1):24*lastDay);
        annualProfile = Project.Result.eOut.Grid_purchases_hourly_kWh.*SO2 + RESULTS.eOut.SO2 + BoilerSO2*RESULTS.eOut.BoilerHeatHour;
        Y2 = annualProfile(1+24*(day1-1):24*lastDay);
        Ylab = 'SO2 Emissions (lb/hr)';
        LegendText = {'Baseline';'With CHP';};
    case 'NOx Emissions'
        baselineProfile = BaselinePurchaseHour.*NOx' + BoilerNOx*BaselineHeatHour;
        Y = baselineProfile(1+24*(day1-1):24*lastDay);
        annualProfile = Project.Result.eOut.Grid_purchases_hourly_kWh.*NOx + RESULTS.eOut.NOx + BoilerNOx*RESULTS.eOut.BoilerHeatHour;
        Y2 = annualProfile(1+24*(day1-1):24*lastDay);
        Ylab = 'NOx Emissions (lb/hr)';
        LegendText = {'Baseline';'With CHP';};
    case 'Grid emission factors'
        Y = [NOx(1+24*(day1-1):24*lastDay) SO2(1+24*(day1-1):24*lastDay)];
        Y2 = CO2(1+24*(day1-1):24*lastDay)./1000;
        Ylab = 'Grid Emission factors (lb/hr or ton/hr)';
        LegendText = {'NOx';'SO2';'CO2';};
    end
plot(X,Y)
hold on
plot(X,Y2,'r');
hold off
ylabel(Ylab)
legend(LegendText)
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



% --- Executes on button press in pushbuttonOptSize.
function pushbuttonOptSize_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonOptSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
MSG = msgbox('Re-sizing CHP System');
ScaleSystemComponents(100,100,100,120,1,2,1,2,2)
updateViewResults(handles,0)
close(MSG)


% --- Executes on selection change in popupmenuControl.
function popupmenuControl_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuControl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Project Model_dir
val=get(hObject,'value');
str=get(hObject,'string');
currentstr=str{val};
Control_dir=fullfile(Model_dir, 'System Library' , 'Control'); %strrep(which('NREL_FCModel.m'),'\main\NREL_FCModel.m','\component library\Control');
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

% --- Executes on selection change in popupmenuSubregion.
function popupmenuSubregion_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuSubregion (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
updateViewResults(handles,1)

% --- Executes during object creation, after setting all properties.
function popupmenuSubregion_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuSubregion (see GCBO)
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
ZipCode = str2double(get(hObject,'String'));
StateNum = get(handles.popupmenuState,'Value');
stateAbrev = {'AL';'AK';'AZ';'AR';'CA';'CO';'CT';'DE';'FL';'GA';'HI';'ID';'IL';'IN';'IA';'KS';'KY';'LA';'ME';'MD';'MA';'MI';'MN';'MS';'MO';
              'MT';'NE';'NV';'NH';'NJ';'NM';'NY';'NC';'ND';'OH';'OK';'OR';'PA';'RI';'SC';'SD';'TN';'TX';'UT';'VT';'VA';'WA';'WV';'WI';'WY';};
State = stateAbrev(StateNum);
ASHRAEZone = ClimateZone( char(State),ZipCode );
EmissionRegion = EmissionZone( char(State),ZipCode );
eGridZones = get(handles.popupmenuSubregion,'string');
regionVal = find(strcmp(EmissionRegion,eGridZones));
set(handles.popupmenuSubregion,'value',regionVal);
Climates = {'1A'; '2A'; '2B'; '3A'; '3B'; '3B-Coast'; '3C'; '4A'; '4B'; '4C'; '5A'; '5B'; '6A'; '6B'; '7'; '8';};
climateVal = find(strcmp(ASHRAEZone,Climates));
set(handles.popupmenuClimate,'value',climateVal);
if climateVal~=get(handles.popupmenuClimate,'value');
    popupmenuBuilding_Callback(hObject, eventdata, handles)
end


% --- Executes during object creation, after setting all properties.
function editZipCode_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editZipCode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
