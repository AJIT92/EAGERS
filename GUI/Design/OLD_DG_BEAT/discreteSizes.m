function varargout = discreteSizes(varargin)
% DISCRETESIZES MATLAB code for discreteSizes.fig
%      DISCRETESIZES, by itself, creates a new DISCRETESIZES or raises the existing
%      singleton*.
%
%      H = DISCRETESIZES returns the handle to a new DISCRETESIZES or the handle to
%      the existing singleton*.
%
%      DISCRETESIZES('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DISCRETESIZES.M with the given input arguments.
%
%      DISCRETESIZES('Property','Value',...) creates a new DISCRETESIZES or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before discreteSizes_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to discreteSizes_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help discreteSizes

% Last Modified by GUIDE v2.5 18-Mar-2014 14:45:04

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @discreteSizes_OpeningFcn, ...
                   'gui_OutputFcn',  @discreteSizes_OutputFcn, ...
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


% --- Executes just before discreteSizes is made visible.
function discreteSizes_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to discreteSizes (see VARARGIN)

% Choose default command line output for discreteSizes
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes discreteSizes wait for user response (see UIRESUME)
uiwait(handles.figure1);



% --- Outputs from this function are returned to the command line.
function varargout = discreteSizes_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure




% --- Executes on button press in pushbuttonYes.
function pushbuttonYes_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonYes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Project
increment = str2double(get(handles.editIncrementSize,'String'));
totalSize = 0;
for i=1:1:length(Project.System.CHP)
    totalSize = totalSize + Project.System.CHP(i).SysSize(1);
end
newSizes = increment*floor(totalSize/increment)/length(Project.System.CHP);

for i=1:1:length(Project.System.CHP)
    Project.System.CHP(i).SysSize(1) = newSizes;
    Project.System.CHP(i).SysSize(2) = newSizes/Project.System.CHP(i).TurnDown;
end
Project.Economic.LifekWh = Project.Economic.LifeYrs*totalSize*8760;
uiresume(handles.figure1);
close(gcf)




% --- Executes on button press in pushbuttonNo.
function pushbuttonNo_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonNo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uiresume(handles.figure1);
close(gcf)




function editIncrementSize_Callback(hObject, eventdata, handles)
% hObject    handle to editIncrementSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editIncrementSize as text
%        str2double(get(hObject,'String')) returns contents of editIncrementSize as a double


% --- Executes during object creation, after setting all properties.
function editIncrementSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editIncrementSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
if(isequal(get(hObject, 'waitstatus'), 'waiting'))
    uiresume(hObject);
else
    delete(hObject);
end
