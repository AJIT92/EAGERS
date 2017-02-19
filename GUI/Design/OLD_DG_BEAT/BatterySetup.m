function varargout = BatterySetup(varargin)
% BATTERYSETUP M-file for BatterySetup.fig
%      BATTERYSETUP, by itself, creates a new BATTERYSETUP or raises the existing
%      singleton*.
%
%      H = BATTERYSETUP returns the handle to a new BATTERYSETUP or the handle to
%      the existing singleton*.
%
%      BATTERYSETUP('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BATTERYSETUP.M with the given input arguments.
%
%      BATTERYSETUP('Property','Value',...) creates a new BATTERYSETUP or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before BatterySetup_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to BatterySetup_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help BatterySetup

% Last Modified by GUIDE v2.5 04-Apr-2013 08:18:02

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @BatterySetup_OpeningFcn, ...
                   'gui_OutputFcn',  @BatterySetup_OutputFcn, ...
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


% --- Executes just before BatterySetup is made visible.
function BatterySetup_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to BatterySetup (see VARARGIN)

% Add a wait dialogue box
figure_Batt = get(0,'CurrentFigure');
MSG_Batt=msgbox('Loading');
figure(figure_Batt)

% Choose default command line output for BatterySetup
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

global Project SYSINDEX BatteryOriginal
BatteryOriginal = Project.System.Battery;
OptimalBatterySizing('editSize');
COMPONENT=Project.System.Battery(SYSINDEX);

set(handles.editName,'string',COMPONENT.Name)
set(handles.textType,'string',COMPONENT.Type)
set(handles.editVersion,'string',COMPONENT.Version)
set(handles.editDescription,'string',COMPONENT.Description)
% set(handles.editPicture,'string',COMPONENT.Picture)
% editPicture_Callback(handles.editPicture,eventdata,handles)
set(handles.editSize,'string',COMPONENT.Size)
set(handles.editSizeMin,'string',COMPONENT.SizeMin)
set(handles.editOptSize,'string',COMPONENT.OptSize)
set(handles.editPeakCharge,'string',COMPONENT.PeakCharge)
set(handles.editPeakDisch,'string',COMPONENT.PeakDisch)
set(handles.editSelfDisch,'string',COMPONENT.SelfDisch)
set(handles.editVoltage,'string',COMPONENT.Voltage)
set(handles.editDOD,'string',COMPONENT.MaxDOD)
set(handles.editChargeResist,'string',COMPONENT.ChargeResist)
set(handles.editDischResist,'string',COMPONENT.DischResist)

if strcmp(COMPONENT.BatteryType, 'leadacid')
    set(handles.uipanelBatteryType,'SelectedObject', handles.radiobuttonLeadAcid)
elseif strcmp(COMPONENT.BatteryType, 'nicad')
    set(handles.uipanelBatteryType,'SelectedObject', handles.radiobuttonNickelCad)
elseif strcmp(COMPONENT.BatteryType, 'nimhyd')
    set(handles.uipanelBatteryType,'SelectedObject', handles.radiobuttonNickleMetHyd)
elseif strcmp(COMPONENT.BatteryType, 'lithium')
    set(handles.uipanelBatteryType,'SelectedObject', handles.radiobuttonLithium)
end

if strcmp(COMPONENT.ChargeMeth, 'CC')
    set(handles.uipanelCharging,'SelectedObject', handles.radiobuttonCC)
elseif strcmp(COMPONENT.ChargeMeth, 'smooth')
    set(handles.uipanelCharging,'SelectedObject', handles.radiobuttonSmooth)
elseif strcmp(COMPONENT.ChargeMeth, 'CV')
    set(handles.uipanelCharging,'SelectedObject', handles.radiobuttonCV)
end

set(handles.uitableBatVol,'Data',COMPONENT.VoltCurve);
axes(findobj(gcf,'tag','axes1'))
plot(Project.System.Battery(SYSINDEX).VoltCurve(:,1),Project.System.Battery(SYSINDEX).VoltCurve(:,2),'b-o');
set(gca,'tag','axes1')
ylim([0 5])
xlabel('% Capacity Discharged')
ylabel('Voltage (V)')

close(MSG_Batt)

% --- Outputs from this function are returned to the command line.
function varargout = BatterySetup_OutputFcn(hObject, eventdata, handles) 
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
Project.System.Battery(SYSINDEX).Name=get(hObject,'string');

% --- Executes during object creation, after setting all properties.
function editName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function editDescription_Callback(hObject, eventdata, handles)
% hObject    handle to editDescription (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Project SYSINDEX
Project.System.Battery(SYSINDEX).Description=get(hObject,'string');

% --- Executes during object creation, after setting all properties.
function editDescription_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editDescription (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function editVersion_Callback(hObject, eventdata, handles)
% hObject    handle to editVersion (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Project SYSINDEX
Project.System.Battery(SYSINDEX).Version=get(hObject,'string');

% --- Executes during object creation, after setting all properties.
function editVersion_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editVersion (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function editPicture_Callback(hObject, eventdata, handles)
% hObject    handle to editPicture (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% global Project SYSINDEX Model_dir
% Project.System.Battery(SYSINDEX).Picture=get(hObject,'string');
% axes(findobj(gcf,'tag','axesImage'));
% image(imread([Model_dir filesep 'graphics' filesep Project.System.Battery(SYSINDEX).Picture]))
% set(gca,'tag','axesImage')
% axis off

% --- Executes during object creation, after setting all properties.
function editPicture_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editPicture (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbuttonPicture.
function pushbuttonPicture_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonPicture (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Project SYSINDEX
[f,~]=uigetfile(Project.System.Battery(SYSINDEX).Picture,'Select Picture File',which(Project.System.Battery(SYSINDEX).Picture));
if f==0; return; end
Project.System.Battery(SYSINDEX).Picture=f;
set(handles.editPicture,'string',Project.System.Battery(SYSINDEX).Picture)
% editPicture_Callback(handles.editPicture,eventdata,handles)

% --- Executes on button press in pushbuttonSaveOnlyToProject.
function pushbuttonSaveOnlyToProject_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonSaveOnlyToProject (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clear global SYSINDEX BatteryOriginal
close(gcf)
uiresume

% --- Executes on button press in pushbuttonCancel.
function pushbuttonCancel_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonCancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Project BatteryOriginal
Project.System.Battery = BatteryOriginal;
clear global SYSINDEX BatteryOriginal
close(gcf)

% --- Executes on button press in pushbuttonSaveAs.
function pushbuttonSaveAs_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonSaveAs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Project SYSINDEX Model_dir
[f,p]=uiputfile('*.mat','Save As Battery Component',fullfile(Model_dir,'System Library','Battery'));
if f==0;return;end
component=Project.System.Battery(SYSINDEX);
save([p f],'component')
clear global SYSINDEX BatteryOriginal
close(gcf)
uiresume

% --- Executes during object creation, after setting all properties.
function uipanelBatteryType_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uipanelBatteryType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



% --- Executes when selected object is changed in uipanelBatteryType.
function uipanelBatteryType_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanelBatteryType 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
global Project SYSINDEX
switch get(eventdata.NewValue,'Tag')
    case 'radiobuttonLeadAcid'
        Project.System.Battery(SYSINDEX).BatteryType = 'leadacid';
    case 'radiobuttonNickelCad'
        Project.System.Battery(SYSINDEX).BatteryType = 'nicad';
    case 'radiobuttonNickleMetHyd'
        Project.System.Battery(SYSINDEX).BatteryType = 'nimhyd';
    case 'radiobuttonLithium'
        Project.System.Battery(SYSINDEX).BatteryType = 'lithium';
end

% --- Executes when selected object is changed in uipanelCharging.
function uipanelCharging_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanelCharging 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
global Project SYSINDEX
switch get(eventdata.NewValue,'Tag')
    case 'radiobuttonCC'
        Project.System.Battery(SYSINDEX).ChargeMeth = 'CC';
    case 'radiobuttonCV'
        Project.System.Battery(SYSINDEX).ChargeMeth = 'CV';
    case 'radiobuttonSmooth'
        Project.System.Battery(SYSINDEX).ChargeMeth = 'smooth';
end

function editSize_Callback(hObject, eventdata, handles)
% hObject    handle to editSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Project SYSINDEX
Project.System.Battery(SYSINDEX).Size=str2double(get(hObject,'string'));
OptimalBatterySizing('editSize');
COMPONENT = Project.System.Battery(SYSINDEX);
set(handles.editSize,'string',COMPONENT.Size)
set(handles.editOptSize,'string',COMPONENT.OptSize)
set(handles.editSizeMin,'string',COMPONENT.SizeMin)


% --- Executes during object creation, after setting all properties.
function editSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function editOptSize_Callback(hObject, eventdata, handles)
% hObject    handle to editOptSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Project SYSINDEX
Project.System.Battery(SYSINDEX).OptSize=str2double(get(hObject,'string'));
OptimalBatterySizing('OptSize');
COMPONENT = Project.System.Battery(SYSINDEX);
set(handles.editSize,'string',COMPONENT.Size)
set(handles.editOptSize,'string',COMPONENT.OptSize)
set(handles.editSizeMin,'string',COMPONENT.SizeMin)

% --- Executes during object creation, after setting all properties.
function editOptSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editOptSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editSizeMin_Callback(hObject, eventdata, handles)
% hObject    handle to editSizeMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Project SYSINDEX
Project.System.Battery(SYSINDEX).SizeMin = str2double(get(hObject,'string'));
OptimalBatterySizing('editSizeMin');
COMPONENT = Project.System.Battery(SYSINDEX);
set(handles.editSize,'string',COMPONENT.Size)
set(handles.editOptSize,'string',COMPONENT.OptSize)
set(handles.editSizeMin,'string',COMPONENT.SizeMin)

% --- Executes during object creation, after setting all properties.
function editSizeMin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editSizeMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function editPeakCharge_Callback(hObject, eventdata, handles)
% hObject    handle to editPeakCharge (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Project SYSINDEX
Project.System.Battery(SYSINDEX).PeakCharge = str2double(get(hObject,'string'));


% --- Executes during object creation, after setting all properties.
function editPeakCharge_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editPeakCharge (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function editPeakDisch_Callback(hObject, eventdata, handles)
% hObject    handle to editPeakDisch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Project SYSINDEX
Project.System.Battery(SYSINDEX).PeakDisch = str2double(get(hObject,'string'));
editSize_Callback(handles.editSize, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function editPeakDisch_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editPeakDisch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function editVoltage_Callback(hObject, eventdata, handles)
% hObject    handle to editVoltage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Project SYSINDEX
Project.System.Battery(SYSINDEX).Voltage = str2double(get(hObject,'string'));
editSize_Callback(handles.editSize, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function editVoltage_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editVoltage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function editDOD_Callback(hObject, eventdata, handles)
% hObject    handle to editDOD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Project SYSINDEX
Project.System.Battery(SYSINDEX).MaxDOD = str2double(get(hObject,'string'));
editSize_Callback(handles.editSize, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function editDOD_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editDOD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function editChargeResist_Callback(hObject, eventdata, handles)
% hObject    handle to editChargeResist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Project SYSINDEX
Project.System.Battery(SYSINDEX).ChargeResist = str2double(get(hObject,'string'));

% --- Executes during object creation, after setting all properties.
function editChargeResist_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editChargeResist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function editDischResist_Callback(hObject, eventdata, handles)
% hObject    handle to editDischResist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Project SYSINDEX
Project.System.Battery(SYSINDEX).DischResist = str2double(get(hObject,'string'));
editSize_Callback(handles.editSize, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function editDischResist_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editDischResist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function editSelfDisch_Callback(hObject, eventdata, handles)
% hObject    handle to editSelfDisch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Project SYSINDEX
Project.System.Battery(SYSINDEX).SelfDisch = str2double(get(hObject,'string'));

% --- Executes during object creation, after setting all properties.
function editSelfDisch_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editSelfDisch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function uitableBatVol_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to uitableEfficiency (see GCBO)
% eventdata  structure with the following fields (see UITABLEEFFICIENCY)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
global Project SYSINDEX
Project.System.Battery(SYSINDEX).VoltCurve = get(handles.uitableBatVol,'Data');
axes(findobj(gcf,'tag','axes1'))
plot(Project.System.Battery(SYSINDEX).VoltCurve(:,1),Project.System.Battery(SYSINDEX).VoltCurve(:,2),'b-o');
set(gca,'tag','axes1')
ylim([0 5])
xlabel('% Capacity Discharged')
ylabel('Voltage (V)')
