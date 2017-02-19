function  varargout = refrigerationChooser(varargin)
% REFRIGERATIONCHOOSER MATLAB code for refrigerationChooser.fig
%      REFRIGERATIONCHOOSER, by itself, creates a new REFRIGERATIONCHOOSER or raises the existing
%      singleton*.
%
%      H = REFRIGERATIONCHOOSER returns the handle to a new REFRIGERATIONCHOOSER or the handle to
%      the existing singleton*.
%
%      REFRIGERATIONCHOOSER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in REFRIGERATIONCHOOSER.M with the given input arguments.
%
%      REFRIGERATIONCHOOSER('Property','Value',...) creates a new REFRIGERATIONCHOOSER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before refrigerationChooser_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to refrigerationChooser_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help refrigerationChooser

% Last Modified by GUIDE v2.5 17-Mar-2014 15:16:05

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @refrigerationChooser_OpeningFcn, ...
                   'gui_OutputFcn',  @refrigerationChooser_OutputFcn, ...
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


% --- Executes just before refrigerationChooser is made visible.
function refrigerationChooser_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to refrigerationChooser (see VARARGIN)

% Choose default command line output for refrigerationChooser
global Model_dir tempBuild name DistCool
handles.output = hObject;
name = char(varargin{1});
buildDir = fullfile(Model_dir,'System Library','Buildings');
load(fullfile(buildDir,name));
tempBuild = component;
DistCool = varargin{2};
plotlist={'Electric Demand'
        'Cooling Demand'};
set(handles.sliderRefrigPerc,'Min',1,'Max',100,'Value',50)
set(handles.editAvgCOP,'string',num2str(2.25));
%note that the original COP's are 3 and 1.5
set(handles.editFreezeCOP,'string',num2str(1.5));
set(handles.editCoolerCOP,'string',num2str(3));
set(handles.editFreezePerc,'string',num2str(50));
set(handles.editCoolerPerc,'string',num2str(50));

loadbuilding(handles)

set(handles.popupmenuAxes,'string',plotlist,'value',1)
set(handles.sliderZoom,'Min',1,'Max',4,'Value',1)
sliderZoom_Callback(handles.sliderZoom, eventdata, handles)


% Update handles structure
guidata(hObject, handles);

% UIWAIT makes refrigerationChooser wait for user response (see UIRESUME)


% --- Outputs from this function are returned to the command line.
function varargout = refrigerationChooser_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global tempBuild
Refrig.ElecLoad = tempBuild.Refrig.ElecLoad;
Refrig.CoolLoad = tempBuild.Refrig.CoolLoad;
varargout{1} = Refrig;

function loadbuilding(handles)
global tempBuild DistCool 
FreezeCOP = str2double(get(handles.editFreezeCOP,'string'));
CoolerCOP = str2double(get(handles.editCoolerCOP,'string'));
FreezePerc = str2double(get(handles.editFreezePerc,'string'));
%note that the original COP's are 3 and 1.5
FreezeLoad = tempBuild.RefrigerationElecLoad*(FreezePerc/100)*1.5/FreezeCOP;
CoolerLoad = tempBuild.RefrigerationElecLoad*(1-FreezePerc/100)*3/CoolerCOP;
tempBuild.Refrig.ElecLoad = (FreezeLoad+CoolerLoad)*get(handles.sliderRefrigPerc,'value')/100;
tempBuild.Refrig.CoolLoad = (FreezeLoad*FreezeCOP+CoolerLoad*CoolerCOP);
tempBuild.DemandC = DistCool*tempBuild.DemandC;
tempBuild.CoolingElectricalLoad = DistCool*tempBuild.CoolingElectricalLoad + tempBuild.Refrig.ElecLoad;

% --- Executes on button press in pushbuttonDone.
function pushbuttonDone_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonDone (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close(gcf)
uiresume


function popupmenuAxes_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuAxes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global tempBuild
Ts = 8760/length(tempBuild.DemandE);
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

val=get(hObject,'value');
str=get(hObject,'string');
currentstr=str{val};
switch currentstr
    case 'Electric Demand'
        PlotType = 'plot';
        Y = tempBuild.DemandE(1+24/Ts*(day1-1):24/Ts*lastDay);
        Y2 = tempBuild.DemandE(1+24/Ts*(day1-1):24/Ts*lastDay) - tempBuild.CoolingElectricalLoad(1+24/Ts*(day1-1):24/Ts*lastDay);
        LegendText = {'Baseline Elecric Demand','Electric Demand with CCHP replacing refrigeration'};  
    case 'Cooling Demand'
        PlotType = 'plot';
        Y = tempBuild.DemandC(1+24/Ts*(day1-1):24/Ts*lastDay);
        Y2 = tempBuild.DemandC(1+24/Ts*(day1-1):24/Ts*lastDay) + tempBuild.Refrig.CoolLoad(1+24/Ts*(day1-1):24/Ts*lastDay);
        LegendText = {'Baseline Cooling Demand','Cooling Demand with CCHP replacing refrigeration'};
end

plot(X,Y);
hold on
plot(X,Y2,'r');
ylabel('Demand (kW)') 
legend(LegendText)
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


function editAvgCOP_Callback(hObject, eventdata, handles)
% hObject    handle to editAvgCOP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
AvgCOP = str2double(get(handles.editAvgCOP,'string'));
set(handles.editAvgCOP,'string',num2str(AvgCOP));
FreezeCOP = str2double(get(handles.editFreezeCOP,'string'));
FreezePerc = str2double(get(handles.editFreezePerc,'string'));
CoolerCOP = str2double(get(handles.editCoolerCOP,'string'));
oldCOP = FreezePerc/100*FreezeCOP+(1-FreezePerc/100)*CoolerCOP;
newFreeze = AvgCOP/oldCOP*FreezeCOP;
newCooler = AvgCOP/oldCOP*CoolerCOP;
set(handles.editFreezeCOP,'string',num2str(newFreeze));
set(handles.editCoolerCOP,'string',num2str(newCooler));
loadbuilding(handles)


% --- Executes during object creation, after setting all properties.
function editAvgCOP_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editAvgCOP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function editCoolerPerc_Callback(hObject, eventdata, handles)
% hObject    handle to editCoolerPerc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
CoolerPerc = str2double(get(handles.CoolerPerc,'string'));
set(handles.editFreezerPerc,'string',num2str(100-CoolerPerc));
FreezeCOP = str2double(get(handles.editFreezeCOP,'string'));
FreezePerc = str2double(get(handles.editFreezePerc,'string'));
CoolerCOP = str2double(get(handles.editCoolerCOP,'string'));
AvgCOP = FreezePerc/100*FreezeCOP+(1-FreezePerc/100)*CoolerCOP;
set(handles.editAvgCOP,'string',num2str(AvgCOP));
loadbuilding(handles)

% --- Executes during object creation, after setting all properties.
function editCoolerPerc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editCoolerPerc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editCoolerCOP_Callback(hObject, eventdata, handles)
% hObject    handle to editCoolerCOP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
CoolerCOP = str2double(get(handles.editCoolerCOP,'string'));
FreezeCOP = str2double(get(handles.editFreezeCOP,'string'));
FreezePerc = str2double(get(handles.editFreezePerc,'string'));
AvgCOP = FreezePerc/100*FreezeCOP+(1-FreezePerc/100)*CoolerCOP;
set(handles.editFreezeCOP,'string',num2str(newFreeze));
set(handles.editCoolerCOP,'string',num2str(newCooler));
set(handles.editAvgCOP,'string',num2str(AvgCOP));
loadbuilding(handles)

% --- Executes during object creation, after setting all properties.
function editCoolerCOP_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editCoolerCOP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function sliderRefrigPerc_Callback(hObject, eventdata, handles)
% hObject    handle to sliderRefrigPerc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
loadbuilding(handles)


% --- Executes during object creation, after setting all properties.
function sliderRefrigPerc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderRefrigPerc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function editFreezePerc_Callback(hObject, eventdata, handles)
% hObject    handle to editFreezePerc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
FreezePerc = str2double(get(handles.editFreezePerc,'string'));
set(handles.editCoolerPerc,'string',num2str(100-FreezePerc));
FreezeCOP = str2double(get(handles.editFreezeCOP,'string'));
CoolerCOP = str2double(get(handles.editCoolerCOP,'string'));
AvgCOP = FreezePerc/100*FreezeCOP+(1-FreezePerc/100)*CoolerCOP;
set(handles.editAvgCOP,'string',num2str(AvgCOP));
loadbuilding(handles)

% --- Executes during object creation, after setting all properties.
function editFreezePerc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editFreezePerc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editFreezeCOP_Callback(hObject, eventdata, handles)
% hObject    handle to editFreezeCOP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
CoolerCOP = str2double(get(handles.editCoolerCOP,'string'));
FreezeCOP = str2double(get(handles.editFreezeCOP,'string'));
FreezePerc = str2double(get(handles.editFreezePerc,'string'));
AvgCOP = FreezePerc/100*FreezeCOP+(1-FreezePerc/100)*CoolerCOP;
set(handles.editFreezeCOP,'string',num2str(newFreeze));
set(handles.editCoolerCOP,'string',num2str(newCooler));
set(handles.editAvgCOP,'string',num2str(AvgCOP));
loadbuilding(handles)

% --- Executes during object creation, after setting all properties.
function editFreezeCOP_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editFreezeCOP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
