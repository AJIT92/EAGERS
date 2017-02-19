function varargout = EditTable(varargin)
% EDITTABLE MATLAB code for EditTable.fig
%      EDITTABLE, by itself, creates a new EDITTABLE or raises the existing
%      singleton*.
%
%      H = EDITTABLE returns the handle to a new EDITTABLE or the handle to
%      the existing singleton*.
%
%      EDITTABLE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EDITTABLE.M with the given input arguments.
%
%      EDITTABLE('Property','Value',...) creates a new EDITTABLE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before EditTable_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to EditTable_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help EditTable

% Last Modified by GUIDE v2.5 09-Oct-2014 11:48:34

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @EditTable_OpeningFcn, ...
                   'gui_OutputFcn',  @EditTable_OutputFcn, ...
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


% --- Executes just before EditTable is made visible.
function EditTable_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to EditTable (see VARARGIN)

% Choose default command line output for EditTable
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

global EditTableHandle TableData
EditTableHandle = handles.figure1;
ColNames = varargin{1};
Text1 = varargin{2};
set(handles.uitable1,'Data',TableData);
set(handles.uitable1,'ColumnName',ColNames);
set(handles.uitable1,'ColumnEditable',true(1,length(TableData(:,1))));
set(handles.text1,'String',Text1);


% --- Outputs from this function are returned to the command line.
function varargout = EditTable_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on button press in pushbuttonRevert.
function pushbuttonRevert_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonRevert (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global TableData
set(handles.uitable1,'Data', TableData)

% --- Executes on button press in pushbuttonFinished.
function pushbuttonFinished_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonFinished (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global EditTableHandle TableData
TableData = get(handles.uitable1,'Data');
close(EditTableHandle)


% --- Executes on button press in pushbuttonAdd.
function pushbuttonAdd_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonAdd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbuttonDelete.
function pushbuttonDelete_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonDelete (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
