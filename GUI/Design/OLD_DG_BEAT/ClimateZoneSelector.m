function varargout = ClimateZoneSelector(varargin)
% CLIMATEZONESELECTOR MATLAB code for ClimateZoneSelector.fig
%      CLIMATEZONESELECTOR, by itself, creates a new CLIMATEZONESELECTOR or raises the existing
%      singleton*.
%
%      H = CLIMATEZONESELECTOR returns the handle to a new CLIMATEZONESELECTOR or the handle to
%      the existing singleton*.
%
%      CLIMATEZONESELECTOR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CLIMATEZONESELECTOR.M with the given input arguments.
%
%      CLIMATEZONESELECTOR('Property','Value',...) creates a new CLIMATEZONESELECTOR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ClimateZoneSelector_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ClimateZoneSelector_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ClimateZoneSelector

% Last Modified by GUIDE v2.5 11-Feb-2014 09:50:49

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ClimateZoneSelector_OpeningFcn, ...
                   'gui_OutputFcn',  @ClimateZoneSelector_OutputFcn, ...
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


% --- Executes just before ClimateZoneSelector is made visible.
function ClimateZoneSelector_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ClimateZoneSelector (see VARARGIN)

% Choose default command line output for ClimateZoneSelector
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
global Project Model_dir
region={
                'AKGD'
                'AKMS'
                'ERCT'
                'FRCC'
                'HIMS'
                'HIOA'
                'MROE'
                'MROW'
                'NYLI'
                'NEWE'
                'NYCW'
                'NYUP'
                'RFCE'
                'RFCM'
                'RFCW'
                'SRMW'
                'SRMV'
                'SRSO'
                'SRTV'
                'SRVC'
                'SPNO'
                'SPSO'
                'CAMX'
                'NWPP'
                'RMPA'
                'AZNM'
              };
set(handles.popupmenuRegion,'string',region,'value',1)

% UIWAIT makes ClimateZoneSelector wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ClimateZoneSelector_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbuttonOK.
function pushbuttonOK_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonOK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Project
popupmenuRegion_Callback(handles.popupmenuRegion,eventdata,handles)
regionVal = get(handles.popupmenuRegion,'value');
Project.Building.Region = regionVal;
close(gcf)

% --- Executes on selection change in popupmenuRegion.
function popupmenuRegion_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuRegion (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenuRegion contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenuRegion


% --- Executes during object creation, after setting all properties.
function popupmenuRegion_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuRegion (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
