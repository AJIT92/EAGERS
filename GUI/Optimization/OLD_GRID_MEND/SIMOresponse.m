function varargout = SIMOresponse(varargin)
% SIMORESPONSE MATLAB code for SIMOresponse.fig
%      SIMORESPONSE, by itself, creates a new SIMORESPONSE or raises the existing
%      singleton*.
%
%      H = SIMORESPONSE returns the handle to a new SIMORESPONSE or the handle to
%      the existing singleton*.
%
%      SIMORESPONSE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SIMORESPONSE.M with the given input arguments.
%
%      SIMORESPONSE('Property','Value',...) creates a new SIMORESPONSE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SIMOresponse_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SIMOresponse_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SIMOresponse

% Last Modified by GUIDE v2.5 23-May-2015 16:20:26

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SIMOresponse_OpeningFcn, ...
                   'gui_OutputFcn',  @SIMOresponse_OutputFcn, ...
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


% --- Executes just before SIMOresponse is made visible.
function SIMOresponse_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SIMOresponse (see VARARGIN)

% Choose default command line output for SIMOresponse
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
global SIMOhandle
SIMOhandle = handles.figure1;
global Generator
A = Generator.VariableStruct.StateSpace.A;
B = Generator.VariableStruct.StateSpace.B;
C = Generator.VariableStruct.StateSpace.C;
D = Generator.VariableStruct.StateSpace.D;
Dt = Generator.VariableStruct.StateSpace.Dt;
SS = ss(A,B,C,D,Dt);
SS2 = d2c(SS);
set(handles.editA,'String',mat2str(SS2(1).A,4));
set(handles.editB,'String',mat2str(SS2(1).B,4));
set(handles.editAdisc,'string',mat2str(A,4))
set(handles.editBdisc,'string',mat2str(B,4))
set(handles.editDt,'String',Dt);
plotResponse(handles,'Discrete')

% --- Outputs from this function are returned to the command line.
function varargout = SIMOresponse_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;

function updateSSMatrices(handles)
global Generator
A = Generator.VariableStruct.StateSpace.A;
B = Generator.VariableStruct.StateSpace.B;
C = Generator.VariableStruct.StateSpace.C;
D = Generator.VariableStruct.StateSpace.D;
Dt = Generator.VariableStruct.StateSpace.Dt;
SS = ss(A,B,C,D,Dt);
SS2 = d2c(SS);
set(handles.editA,'String',mat2str(SS2(1).A,4));
set(handles.editB,'String',mat2str(SS2(1).B,4));
set(handles.editAdisc,'string',mat2str(A,4))
set(handles.editBdisc,'string',mat2str(B,4))
set(handles.editDt,'String',Dt);


function editA_Callback(hObject, eventdata, handles)
plotResponse(handles,'Continuous')

function editB_Callback(hObject, eventdata, handles)
plotResponse(handles,'Continuous')

function editAdisc_Callback(hObject, eventdata, handles)
plotResponse(handles,'Discrete')

function editBdisc_Callback(hObject, eventdata, handles)
plotResponse(handles,'Discrete')


function plotResponse(handles,CorD)
global Generator
Dt = str2num(get(handles.editDt ,'String'));
if strcmp(CorD,'Continuous')
    A = str2num(get(handles.editA,'String'));
    B = str2num(get(handles.editB ,'String'));
    C = [1 0 0 0; 0 0 1 0;];
    D = [0;0;];
    SS = ss(A,B,C,D);
    SS = c2d(SS,1);
    A = SS(1).A;
    B = SS(1).B;
else
    A = str2num(get(handles.editAdisc,'String'));
    B = str2num(get(handles.editBdisc ,'String'));
    C = [1 0 0 0; 0 0 1 0;];
    D = [0;0;];
    SS = ss(A,B,C,D,Dt);
end
Generator.VariableStruct.StateSpace.A = A;
Generator.VariableStruct.StateSpace.B = B;
Generator.VariableStruct.StateSpace.C = C;
Generator.VariableStruct.StateSpace.D = D;
Generator.VariableStruct.StateSpace.Dt = Dt;
x0 = [];
Names ={};
if nnz(Generator.Output.Electricity)>0
    Names(end+1) = {'Electricity'};
    x0(end+1) = Generator.VariableStruct.Startup.Electricity(end);
end
if nnz(Generator.Output.Cooling)>0
    Names(end+1) = {'Cooling'};
    x0(end+1) = Generator.VariableStruct.Startup.Cooling(end);
end
if nnz(Generator.Output.Steam)>0
    Names(end+1) = {'Steam'};
    x0(end+1) = Generator.VariableStruct.Startup.Steam(end);
end
if nnz(Generator.Output.Heat)>0
    Names(end+1) = {'Heat'};
    x0(end+1) = Generator.VariableStruct.Startup.Heat(end);
end
nS = round(3600/Dt)+1;
t = linspace(0, Dt*(nS-1),nS);
u = Generator.Size*linspace(1,1,nS);
[n,n2] = size(C);
X0 = zeros(n2,1);
for i = 1:1:n
    X0(find(C(i,:),1,'first'))=x0(i);
end
[y,t] = lsim(SS,u,t,X0);
axes(handles.axes1);
plot(t/60,y);
xlabel('Time (min)')
legend(Names)
updateSSMatrices(handles)

% --- Executes on button press in pushbuttonClose.
function pushbuttonClose_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonClose (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global SIMOhandle
close(SIMOhandle)

% --- Executes during object creation, after setting all properties.
function editA_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --- Executes during object creation, after setting all properties.
function editB_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --- Executes during object creation, after setting all properties.
function editAdisc_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --- Executes during object creation, after setting all properties.
function editBdisc_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editDt_Callback(hObject, eventdata, handles)
% hObject    handle to editDt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
A = str2num(get(handles.editAdisc,'String'));
B = str2num(get(handles.editBdisc ,'String'));
C = [1 0 0 0; 0 0 1 0;];
D = [0;0;];
Dt = str2num(get(handles.editDt ,'String'));
plotResponse(handles,'Discrete')

% --- Executes during object creation, after setting all properties.
function editDt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editDt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
