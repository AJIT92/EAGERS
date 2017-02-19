function varargout = SetupBatteryStorage(varargin)
% SETUPBATTERYSTORAGE MATLAB code for SetupBatteryStorage.fig
%      SETUPBATTERYSTORAGE, by itself, creates a new SETUPBATTERYSTORAGE or raises the existing
%      singleton*.
%
%      H = SETUPBATTERYSTORAGE returns the handle to a new SETUPBATTERYSTORAGE or the handle to
%      the existing singleton*.
%
%      SETUPBATTERYSTORAGE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SETUPBATTERYSTORAGE.M with the given input arguments.
%
%      SETUPBATTERYSTORAGE('Property','Value',...) creates a new SETUPBATTERYSTORAGE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SetupBatteryStorage_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SetupBatteryStorage_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SetupBatteryStorage

% Last Modified by GUIDE v2.5 14-May-2015 10:17:31

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SetupBatteryStorage_OpeningFcn, ...
                   'gui_OutputFcn',  @SetupBatteryStorage_OutputFcn, ...
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


% --- Executes just before SetupBatteryStorage is made visible.
function SetupBatteryStorage_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SetupBatteryStorage (see VARARGIN)

% Choose default command line output for SetupBatteryStorage
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

global Plant SYSINDEX figHandle Generator
figHandle = handles.figure1;
Generator = Plant.Generator(SYSINDEX);
BAT=Plant.Generator(SYSINDEX).VariableStruct;

set(handles.editName,'string',Generator.Name)
set(handles.editSize,'string',Generator.Size)

set(handles.editPeakCharge,'string',BAT.PeakCharge)
set(handles.editPeakDisch,'string',BAT.PeakDisch)
set(handles.editSelfDisch,'string',BAT.SelfDischarge*(31*24*100));%convert % per hour to %/month)
set(handles.editVoltage,'string',BAT.Voltage)
set(handles.editDOD,'string',BAT.MaxDOD)
set(handles.editChargeResist,'string',BAT.ChargeResist)
set(handles.editDischResist,'string',BAT.DischResist)
set(handles.uitableBatVol,'Data',BAT.VoltCurve);
axes(findobj(gcf,'tag','axes1'))
plot(BAT.VoltCurve(:,1),BAT.VoltCurve(:,2),'b-o');
set(gca,'tag','axes1')
ylim([0 5])
xlabel('% Capacity Discharged')
ylabel('Voltage (V)')


% --- Outputs from this function are returned to the command line.
function varargout = SetupBatteryStorage_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function editSize_Callback(hObject, eventdata, handles)
% hObject    handle to editSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Generator
Generator.Size=str2double(get(hObject,'string'));

% --- Executes during object creation, after setting all properties.
function editSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function editName_Callback(hObject, eventdata, handles)
% hObject    handle to editName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Generator
Generator.Name=get(hObject,'string');


% --- Executes during object creation, after setting all properties.
function editName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbuttonSaveOnlyToPlant.
function pushbuttonSaveOnlyToPlant_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonSaveOnlyToPlant (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Plant SYSINDEX Generator
Plant.Generator(SYSINDEX) = Generator;
clear global Generator
close(gcf)


% --- Executes on button press in pushbuttonCancel.
function pushbuttonCancel_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonCancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close(gcf)


% --- Executes on button press in pushbuttonSaveAs.
function pushbuttonSaveAs_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonSaveAs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Generator Model_dir
if isfield(Generator,'OpMatA')%if the model has already been run
    component = rmfield(Generator,'OpMatA');%remove OpMatA when saving the new component
else
    component=Generator;    
end
savedir=fullfile(Model_dir,'System Library','Electric Storage',strcat(Generator.Name,'.mat'));
[f,p]=uiputfile(savedir,'Save As Energy Storage Component');
if f==0;return;end

save([p f],'component')
pushbuttonSaveOnlyToPlant_Callback(hObject, eventdata, handles)


function editPeakCharge_Callback(hObject, eventdata, handles)
% hObject    handle to editPeakCharge (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Generator
Generator.VariableStruct.PeakCharge = str2double(get(hObject,'string'));

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
global Generator
Generator.VariableStruct.PeakDisch = str2double(get(hObject,'string'));

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
global Generator
Generator.VariableStruct.Voltage = str2double(get(hObject,'string'));
 
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
global Generator
Generator.VariableStruct.MaxDOD = str2double(get(hObject,'string'));
 

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
global Generator
Generator.VariableStruct.ChargeResist = str2double(get(hObject,'string'));

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
global Generator
Generator.VariableStruct.DischResist = str2double(get(hObject,'string'));
 
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
global Generator
Generator.VariableStruct.SelfDischarge = str2double(get(hObject,'string'))/(31*24*100);%convert %/month to fraction of SOC per hour

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
global Generator
Generator.VariableStruct.VoltCurve = get(handles.uitableBatVol,'Data');
axes(findobj(gcf,'tag','axes1'))
plot(Generator.VariableStruct.VoltCurve(:,1),Generator.VariableStruct.VoltCurve(:,2),'b-o');
set(gca,'tag','axes1')
ylim([0 5])
xlabel('% Capacity Discharged')
ylabel('Voltage (V)')


% --- Executes on button press in pushbuttonComPorts.
function pushbuttonComPorts_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonComPorts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global EditCommHandle
CommunicationPorts();
waitfor(EditCommHandle)