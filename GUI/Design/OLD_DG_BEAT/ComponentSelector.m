function varargout = ComponentSelector(varargin)
% COMPONENTSELECTOR MATLAB code for ComponentSelector.fig
%      COMPONENTSELECTOR, by itself, creates a new COMPONENTSELECTOR or raises the existing
%      singleton*.
%
%      H = COMPONENTSELECTOR returns the handle to a new COMPONENTSELECTOR or the handle to
%      the existing singleton*.
%
%      COMPONENTSELECTOR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in COMPONENTSELECTOR.M with the given input arguments.
%
%      COMPONENTSELECTOR('Property','Value',...) creates a new COMPONENTSELECTOR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ComponentSelector_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ComponentSelector_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ComponentSelector

% Last Modified by GUIDE v2.5 19-Jul-2013 14:19:08

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ComponentSelector_OpeningFcn, ...
                   'gui_OutputFcn',  @ComponentSelector_OutputFcn, ...
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
% End initialization code - DO NOT EDIT THIS LINE


% --- Executes just before ComponentSelector is made visible.
function ComponentSelector_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ComponentSelector (see VARARGIN)

% Choose default command line output for ComponentSelector
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ComponentSelector wait for user response (see UIRESUME)
% uiwait(handles.figure1);

list={'CHP'
    'Chiller'
    'Battery'
    'ThermalStorage'
    'Solar'
    'Wind'};
set(handles.popupmenuType,'string',list,'value',1)
popupmenuType_Callback(handles.popupmenuType,eventdata,handles)


% --- Outputs from this function are returned to the command line.
function varargout = ComponentSelector_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in popupmenuType.
function popupmenuType_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global ComponentSelectorPick Model_dir
v=get(hObject,'value');
s=get(hObject,'string');
ComponentSelectorPick.Type=s{v};
compdir=fullfile(Model_dir, 'System Library', ComponentSelectorPick.Type); %fullfile(strrep(which('NREL_FCModel.m'),'\main\NREL_FCModel.m','\component library'), ComponentSelectorPick.Type);
files=dir(fullfile(compdir,'*.mat'));
list2=strrep({files.name},'.mat','');
if strcmp(ComponentSelectorPick.Type,'Solar')
    n = 0;
    for i = 1:1:length(list2)
        if (strcmp(list2(i),'DirectNormal')+strcmp(list2(i),'GlobalHorizontal')+strcmp(list2(i),'Azimuth')+strcmp(list2(i),'Zenith')) == 0
            n = n+1;
            list3(n) = list2(i);
        end
    end
    list2 = list3;
elseif strcmp(ComponentSelectorPick.Type,'Wind')
    n = 0;
    for i = 1:1:length(list2)
        if (strcmp(list2(i),'WindPow')+strcmp(list2(i),'WindSpeed')) == 0
            n = n+1;
            list3(n) = list2(i);
        end
    end
    list2 = list3;
end
ComponentSelectorPick.Selection=list2{1};
set(handles.popupmenuSelection,'string',list2,'value',1)


% --- Executes during object creation, after setting all properties.
function popupmenuType_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenuSelection.
function popupmenuSelection_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuSelection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global ComponentSelectorPick
s=get(hObject,'string');
v=get(hObject,'value');
ComponentSelectorPick.Selection=s{v};

% --- Executes during object creation, after setting all properties.
function popupmenuSelection_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuSelection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbuttonOK.
function pushbuttonOK_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonOK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Project ComponentSelectorPick Model_dir
% filename=fullfile(strrep(which('NREL_FCModel.m'),'\main\NREL_FCModel.m','\component library'), ComponentSelectorPick.Type,ComponentSelectorPick.Selection);
filename = fullfile(Model_dir, 'System Library',ComponentSelectorPick.Type,ComponentSelectorPick.Selection);
load(filename)
if strcmp(ComponentSelectorPick.Type,'Solar') || strcmp(ComponentSelectorPick.Type,'Wind')
    if ~isfield(Project,'Renewable')
        Project.Renewable.Name = 'Renewables';
    end
    if isfield(Project.Renewable,component.Type)
        Project.Renewable.(component.Type)(end+1)=component;
    else Project.Renewable.(component.Type)=component;
    end
elseif isfield(Project.System,component.Type)
    Project.System.(component.Type)(end+1)=component;
else Project.System.(component.Type)=component;
end
close(gcf)
clear global ComponentSelectorPick

% --- Executes on button press in pushbuttonCancel.
function pushbuttonCancel_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonCancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clear global ComponentSelectorPick
close(gcf)
