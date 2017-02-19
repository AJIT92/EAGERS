function varargout = WindSetup(varargin)
% WINDSETUP MATLAB code for WindSetup.fig
%      WINDSETUP, by itself, creates a new WINDSETUP or raises the existing
%      singleton*.
%
%      H = WINDSETUP returns the handle to a new WINDSETUP or the handle to
%      the existing singleton*.
%
%      WINDSETUP('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in WINDSETUP.M with the given input arguments.
%
%      WINDSETUP('Property','Value',...) creates a new WINDSETUP or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before WindSetup_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to WindSetup_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help WindSetup

% Last Modified by GUIDE v2.5 23-Apr-2013 08:57:12

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @WindSetup_OpeningFcn, ...
                   'gui_OutputFcn',  @WindSetup_OutputFcn, ...
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


% --- Executes just before WindSetup is made visible.
function WindSetup_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Add a wait dialogue box
figure_Wind = get(0,'CurrentFigure');
MSG_Wind=msgbox('Loading');
figure(figure_Wind)

% varargin   command line arguments to WindSetup (see VARARGIN)
global Project SYSINDEX  WindOriginal
WindOriginal = Project.Renewable.Wind(SYSINDEX);
% Choose default command line output for WindSetup
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
wind = Project.Renewable.Wind(SYSINDEX);
stateName = {'Alabama';'Alaska';'Arizona';'Arkansas';'California';'Colorado';'Connecticut';'Delaware';'Florida';'Georgia';'Hawaii';'Idaho';'Illinois';'Indiana';'Iowa';'Kansas';
             'Kentucky';'Louisiana';'Maine';'Maryland';'Massachusetts';'Michigan';'Minnesota';'Mississippi';'Missouri';'Montana';'Nebraska';'Nevada';'New Hampshire';'New Jersey';
             'New Mexico';'New York';'North Carolina';'North Dakota';'Ohio';'Oklahoma';'Oregon';'Pennsylvania';'Rhode Island';'South Carolina';'South Dakota';'Tennessee';'Texas';
             'Utah';'Vermont';'Virginia';'Washington';'West Virginia';'Wisconsin';'Wyoming';};
stateNum = find(strcmp(wind.State,stateName));

set(handles.popupmenuState,'string',stateName,'value',stateNum)        
set(handles.editName,'string',wind.Name)
set(handles.editSize,'string',wind.Size)
set(handles.editEfficiency,'string',wind.Eff)
set(handles.editDiameter,'string',wind.Diam)
set(handles.editHeight,'string',wind.Height)
set(handles.editElevation,'string',wind.Elev)
set(handles.editCutIn,'string',wind.CutIn)
set(handles.editRated,'string',wind.Rated)
set(handles.editShutDown,'string',wind.ShutDown)

wind = Project.Renewable.Wind;
CHPsize = 0;
for i = 1:1:length(Project.System.CHP)
    CHPsize = CHPsize + Project.System.CHP(i).SysSize(1);
end
AnnualGen = 8760*CHPsize;
AnnualDem = sum(Project.Building.DemandE)*8760/length(Project.Building.DemandE);
AnnualPower = zeros(length(wind),1);
for i = 1:1:length(wind)
    AnnualPower(i) = annualPower(Project.Renewable.Wind(i));
    wind(i).GenFrac = AnnualPower(i)/AnnualGen*100;
    wind(i).DemFrac = AnnualPower(i)/AnnualDem*100;
end
set(handles.editSizeAgen,'string',wind(SYSINDEX).GenFrac)
set(handles.editSizeAdem,'string',wind(SYSINDEX).DemFrac)
Project.Renewable.Wind = wind;

close(MSG_Wind)

% --- Outputs from this function are returned to the command line.
function varargout = WindSetup_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function editName_Callback(hObject, eventdata, handles)
% hObject    handle to editName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Project SYSINDEX
Project.Renewable.Wind(SYSINDEX).Name=get(hObject,'string');

% --- Executes during object creation, after setting all properties.
function editName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editName (see GCBO)
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
global Project SYSINDEX Model_dir
val=get(hObject,'value');
str=get(hObject,'string');
Project.Renewable.Wind(SYSINDEX).State=str{val};
WindDir = fullfile (Model_dir, 'System Library' , 'Wind'); % strrep(which('NREL_FCModel.m'),'\main\NREL_FCModel.m','\component library\Wind');
load(strcat(WindDir,[filesep 'windData' filesep 'WindSpeed.mat']));
load(strcat(WindDir,[filesep 'windData' filesep 'WindPow.mat']));
Project.Renewable.Wind(SYSINDEX).Wind = WindSpeed(:,val);
CHPsize = 0;
for i = 1:1:length(Project.System.CHP)
    CHPsize = CHPsize + Project.System.CHP(i).SysSize(1);
end
AnnualGen = 8760*CHPsize;
AnnualDem = sum(Project.Building.DemandE)*8760/length(Project.Building.DemandE);
AnnualPower = annualPower(Project.Renewable.Wind(SYSINDEX));
wind(SYSINDEX).GenFrac = AnnualPower/AnnualGen*100;
wind(SYSINDEX).DemFrac = AnnualPower/AnnualDem*100;
set(handles.editSizeAgen,'string',wind(SYSINDEX).GenFrac)
set(handles.editSizeAdem,'string',wind(SYSINDEX).DemFrac)

% --- Executes during object creation, after setting all properties.
function popupmenuState_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuState (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editSize_Callback(hObject, eventdata, handles)
% hObject    handle to editSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Project SYSINDEX
wind = Project.Renewable.Wind;
newSize = str2double(get(hObject,'string'));
wind(SYSINDEX).Size = newSize;
rho = 1.1798-1.3793e-4*wind(SYSINDEX).Elev+5.667e-9*wind(SYSINDEX).Elev^2;
Eff = (newSize*1000)/wind(SYSINDEX).Rated.^3*(8/(3.1416*rho*wind(SYSINDEX).Diam^2));
if Eff>.45 || Eff<.2
    diam = (newSize*1000)/wind(SYSINDEX).Rated.^3*(8/(3.1416*rho*wind(SYSINDEX).Eff));
    wind(SYSINDEX).Diam = diam;
    set(handles.editDiameter,'string',wind(SYSINDEX).Diam)
end
set(handles.editEfficiency,'string',wind(SYSINDEX).Eff)
CHPsize = 0;
for i = 1:1:length(Project.System.CHP)
    CHPsize = CHPsize + Project.System.CHP(i).SysSize(1);
end
AnnualGen = 8760*CHPsize;
AnnualDem = sum(Project.Building.DemandE)*8760/length(Project.Building.DemandE);
AnnualPower = annualPower(Project.Renewable.Wind(SYSINDEX));
wind(SYSINDEX).GenFrac = AnnualPower/AnnualGen*100;
wind(SYSINDEX).DemFrac = AnnualPower/AnnualDem*100;
set(handles.editSizeAgen,'string',wind(SYSINDEX).GenFrac)
set(handles.editSizeAdem,'string',wind(SYSINDEX).DemFrac)
Project.Renewable.Wind = wind;




% --- Executes during object creation, after setting all properties.
function editSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editEfficiency_Callback(hObject, eventdata, handles)
% hObject    handle to editEfficiency (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Project SYSINDEX
wind = Project.Renewable.Wind;
newEff = str2double(get(hObject,'string'));
wind(SYSINDEX).Eff = newEff;
rho = 1.1798-1.3793e-4*wind(SYSINDEX).Elev+5.667e-9*wind(SYSINDEX).Elev^2;
ratedPow = wind(SYSINDEX).Eff*.5*rho*(3.1416*(wind(SYSINDEX).Diam)^2/4)*wind(SYSINDEX).Rated.^3*(1/1000);
wind(SYSINDEX).Size = ratedPow;
set(handles.editSize,'string',wind(SYSINDEX).Size)

CHPsize = 0;
for i = 1:1:length(Project.System.CHP)
    CHPsize = CHPsize + Project.System.CHP(i).SysSize(1);
end
AnnualGen = 8760*CHPsize;
AnnualDem = sum(Project.Building.DemandE)*8760/length(Project.Building.DemandE);
AnnualPower = annualPower(Project.Renewable.Wind(SYSINDEX));
wind(SYSINDEX).GenFrac = AnnualPower/AnnualGen*100;
wind(SYSINDEX).DemFrac = AnnualPower/AnnualDem*100;
set(handles.editSizeAgen,'string',wind(SYSINDEX).GenFrac)
set(handles.editSizeAdem,'string',wind(SYSINDEX).DemFrac)
Project.Renewable.Wind = wind;

% --- Executes during object creation, after setting all properties.
function editEfficiency_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editEfficiency (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editDiameter_Callback(hObject, eventdata, handles)
% hObject    handle to editDiameter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Project SYSINDEX
wind = Project.Renewable.Wind;
newDiam = str2double(get(hObject,'string'));
wind(SYSINDEX).Diam = newDiam;
rho = 1.1798-1.3793e-4*wind(SYSINDEX).Elev+5.667e-9*wind(SYSINDEX).Elev^2;
ratedPow = wind(SYSINDEX).Eff*.5*rho*(3.1416*(wind(SYSINDEX).Diam)^2/4)*wind(SYSINDEX).Rated.^3*(1/1000);
wind(SYSINDEX).Size = ratedPow;
set(handles.editSize,'string',wind(SYSINDEX).Size)
CHPsize = 0;
for i = 1:1:length(Project.System.CHP)
    CHPsize = CHPsize + Project.System.CHP(i).SysSize(1);
end
AnnualGen = 8760*CHPsize;
AnnualDem = sum(Project.Building.DemandE)*8760/length(Project.Building.DemandE);
AnnualPower = annualPower(Project.Renewable.Wind(SYSINDEX));
wind(SYSINDEX).GenFrac = AnnualPower/AnnualGen*100;
wind(SYSINDEX).DemFrac = AnnualPower/AnnualDem*100;
set(handles.editSizeAgen,'string',wind(SYSINDEX).GenFrac)
set(handles.editSizeAdem,'string',wind(SYSINDEX).DemFrac)
Project.Renewable.Wind = wind;

% --- Executes during object creation, after setting all properties.
function editDiameter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editDiameter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editElevation_Callback(hObject, eventdata, handles)
% hObject    handle to editElevation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Project SYSINDEX
Project.Renewable.Wind(SYSINDEX).Elev=str2double(get(hObject,'string'));
set(handles.editElevation,'string',Project.Renewable.Wind(SYSINDEX).Elev)

% --- Executes during object creation, after setting all properties.
function editElevation_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editElevation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editHeight_Callback(hObject, eventdata, handles)
% hObject    handle to editHeight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Project SYSINDEX
Project.Renewable.Wind(SYSINDEX).Height=str2double(get(hObject,'string'));
set(handles.editHeight,'string',Project.Renewable.Wind(SYSINDEX).Height)

% --- Executes during object creation, after setting all properties.
function editHeight_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editHeight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editRated_Callback(hObject, eventdata, handles)
% hObject    handle to editRated (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Project SYSINDEX
Project.Renewable.Wind(SYSINDEX).Rated=str2double(get(hObject,'string'));
set(handles.editRated,'string',Project.Renewable.Wind(SYSINDEX).Rated)

% --- Executes during object creation, after setting all properties.
function editRated_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editRated (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editCutIn_Callback(hObject, eventdata, handles)
% hObject    handle to editCutIn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Project SYSINDEX
Project.Renewable.Wind(SYSINDEX).CutIn=str2double(get(hObject,'string'));
set(handles.editCutIn,'string',Project.Renewable.Wind(SYSINDEX).CutIn)

% --- Executes during object creation, after setting all properties.
function editCutIn_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editCutIn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editShutDown_Callback(hObject, eventdata, handles)
% hObject    handle to editShutDown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Project SYSINDEX
Project.Renewable.Wind(SYSINDEX).ShutDown=str2double(get(hObject,'string'));
set(handles.editShutDown,'string',Project.Renewable.Wind(SYSINDEX).ShutDown)

% --- Executes during object creation, after setting all properties.
function editShutDown_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editShutDown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pushbuttonSaveOnlyToProject.
function pushbuttonSaveOnlyToProject_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonSaveOnlyToProject (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clear global SYSINDEX WindOriginal
close(gcf)
uiresume

% --- Executes on button press in pushbuttonCancel.
function pushbuttonCancel_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonCancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Project SYSINDEX WindOriginal
Project.Renewable.Wind(SYSINDEX) = WindOriginal;
clear global SYSINDEX WindOriginal
close(gcf)


% --- Executes on button press in pushbuttonSaveAs.
function pushbuttonSaveAs_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonSaveAs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Project SYSINDEX  Model_dir
[f,p]=uiputfile('*.mat','Save As Wind Component',fullfile(Model_dir,'System Library','Wind'));
if f==0;return;end
component=Project.Renewable.Wind(SYSINDEX);
save([p f],'component')
clear global SYSINDEX WindOriginal 
close(gcf)
uiresume


function editSizeAgen_Callback(hObject, eventdata, handles)
% hObject    handle to editSizeAgen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Project SYSINDEX
SizeAgen = str2double(get(hObject,'String'));
WindSizing(SYSINDEX,'generation','no',SizeAgen)
set(handles.editSize,'string',Project.Renewable.Wind(SYSINDEX).Size)
set(handles.editDiameter,'string',Project.Renewable.Wind(SYSINDEX).Diam)
set(handles.editSizeAgen,'string',Project.Renewable.Wind(SYSINDEX).GenFrac)
set(handles.editSizeAdem,'string',Project.Renewable.Wind(SYSINDEX).DemFrac)

% --- Executes during object creation, after setting all properties.
function editSizeAgen_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editSizeAgen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editSizeAdem_Callback(hObject, eventdata, handles)
% hObject    handle to editSizeAdem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Project SYSINDEX
SizeAdem = str2double(get(hObject,'String'));
WindSizing(SYSINDEX,'demand','no',SizeAdem)
set(handles.editSize,'string',Project.Renewable.Wind(SYSINDEX).Size)
set(handles.editDiameter,'string',Project.Renewable.Wind(SYSINDEX).Diam)
set(handles.editSizeAgen,'string',Project.Renewable.Wind(SYSINDEX).GenFrac)
set(handles.editSizeAdem,'string',Project.Renewable.Wind(SYSINDEX).DemFrac)

% --- Executes during object creation, after setting all properties.
function editSizeAdem_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editSizeAdem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function Power = annualPower(wind)
rho = 1.1798-1.3793e-4*wind.Elev+5.667e-9*wind.Elev^2;
CorrectedWind = (wind.Height/80)^.2*wind.Wind;
Wind = CorrectedWind.*(wind.Wind>wind.CutIn).*(wind.Wind<wind.ShutDown);
Wind = min(Wind,wind.Rated);
Power = 8760*wind.Eff*.5*rho*(3.1416*(wind.Diam)^2/4)*sum(Wind.^3)/length(Wind)*(1/1000);
