function varargout = CommunicationPorts(varargin)
% COMMUNICATIONPORTS MATLAB code for CommunicationPorts.fig
%      COMMUNICATIONPORTS, by itself, creates a new COMMUNICATIONPORTS or raises the existing
%      singleton*.
%
%      H = COMMUNICATIONPORTS returns the handle to a new COMMUNICATIONPORTS or the handle to
%      the existing singleton*.
%
%      COMMUNICATIONPORTS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in COMMUNICATIONPORTS.M with the given input arguments.
%
%      COMMUNICATIONPORTS('Property','Value',...) creates a new COMMUNICATIONPORTS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before CommunicationPorts_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to CommunicationPorts_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help CommunicationPorts

% Last Modified by GUIDE v2.5 14-May-2015 10:25:33

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @CommunicationPorts_OpeningFcn, ...
                   'gui_OutputFcn',  @CommunicationPorts_OutputFcn, ...
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


% --- Executes just before CommunicationPorts is made visible.
function CommunicationPorts_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to CommunicationPorts (see VARARGIN)

% Choose default command line output for CommunicationPorts
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
global Generator Comm Measure EditCommHandle
EditCommHandle = handles.figure1;
if isfield(Generator.VariableStruct,'Comm')
    Comm = Generator.VariableStruct.Comm;
else Comm.OnOff = 0;
    Comm.Set = 0;
end
if isfield(Generator.VariableStruct,'Measure')
    Measure = Generator.VariableStruct.Measure;
else Measure.OnOff = 0;
    Measure.Input = 0;
    Measure.Electric = 0;
    Measure.Thermal = 0;
end
set(handles.editCommandOnOff,'string',num2str(Comm.OnOff))
set(handles.editCommandSet,'string',num2str(Comm.Set))
set(handles.editMeasureOnOff,'string',num2str(Measure.OnOff))
set(handles.editMeasureInput,'string',num2str(Measure.Input))
set(handles.editMeasureElectric,'string',num2str(Measure.Electric))
set(handles.editMeasureThermal,'string',num2str(Measure.Thermal))


% --- Outputs from this function are returned to the command line.
function varargout = CommunicationPorts_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function editCommandOnOff_Callback(hObject, eventdata, handles)
% hObject    handle to editCommandOnOff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Comm
Comm.OnOff = str2double(get(hObject,'String'));

% --- Executes during object creation, after setting all properties.
function editCommandOnOff_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editCommandOnOff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editCommandSet_Callback(hObject, eventdata, handles)
% hObject    handle to editCommandSet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Comm
Comm.Set = str2double(get(hObject,'String'));

% --- Executes during object creation, after setting all properties.
function editCommandSet_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editCommandSet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editMeasureThermal_Callback(hObject, eventdata, handles)
% hObject    handle to editMeasureThermal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Measure
Measure.Thermal = str2double(get(hObject,'String'));


% --- Executes during object creation, after setting all properties.
function editMeasureThermal_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editMeasureThermal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function editMeasureElectric_Callback(hObject, eventdata, handles)
% hObject    handle to editMeasureElectric (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Measure
Measure.Electric = str2double(get(hObject,'String'));

% --- Executes during object creation, after setting all properties.
function editMeasureElectric_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editMeasureElectric (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editMeasureInput_Callback(hObject, eventdata, handles)
% hObject    handle to editMeasureInput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Measure
Measure.Input = str2double(get(hObject,'String'));

% --- Executes during object creation, after setting all properties.
function editMeasureInput_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editMeasureInput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editMeasureOnOff_Callback(hObject, eventdata, handles)
% hObject    handle to editMeasureOnOff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Measure
Measure.OnOff = str2double(get(hObject,'String'));

% --- Executes during object creation, after setting all properties.
function editMeasureOnOff_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editMeasureOnOff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pushbuttonFinish.
function pushbuttonFinish_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonFinish (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global EditCommHandle Generator Comm Measure
Generator.VariableStruct.Comm = Comm;
Generator.VariableStruct.Measure = Measure;
close(EditCommHandle)
