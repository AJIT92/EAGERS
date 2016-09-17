function varargout = OptimOptions(varargin)
% OPTIMOPTIONS MATLAB code for OptimOptions.fig
%      OPTIMOPTIONS, by itself, creates a new OPTIMOPTIONS or raises the existing
%      singleton*.
%
%      H = OPTIMOPTIONS returns the handle to a new OPTIMOPTIONS or the handle to
%      the existing singleton*.
%
%      OPTIMOPTIONS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in OPTIMOPTIONS.M with the given input arguments.
%
%      OPTIMOPTIONS('Property','Value',...) creates a new OPTIMOPTIONS or raises
%      the existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before OptimOptions_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to OptimOptions_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help OptimOptions

% Last Modified by GUIDE v2.5 22-Aug-2016 14:16:34

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @OptimOptions_OpeningFcn, ...
                   'gui_OutputFcn',  @OptimOptions_OutputFcn, ...
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

% --- Executes just before OptimOptions is made visible.
function OptimOptions_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to OptimOptions (see VARARGIN)

% Choose default command line output for OptimOptions
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

initialize_gui(hObject, handles, false);

% UIWAIT makes OptimOptions wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --------------------------------------------------------------------
function initialize_gui(fig_handle, handles, isreset)
global Plant 
set(handles.Interval,'string',Plant.optimoptions.Interval);
set(handles.Horizon, 'string', Plant.optimoptions.Horizon);
set(handles.Resolution, 'string', Plant.optimoptions.Resolution);
set(handles.Topt, 'string', Plant.optimoptions.Topt);
set(handles.Tmpc, 'string', Plant.optimoptions.Tmpc);
set(handles.nsSmooth, 'string', Plant.optimoptions.nsSmooth);
set(handles.scaletime, 'string', Plant.optimoptions.scaletime);
set(handles.fastsimulation, 'value', Plant.optimoptions.fastsimulation);
set(handles.slowsimulation, 'value', ~Plant.optimoptions.fastsimulation);
set(handles.constant, 'value', strcmp(Plant.optimoptions.tspacing,'constant'));
set(handles.linear, 'value', strcmp(Plant.optimoptions.tspacing, 'linear'));
set(handles.logarithm, 'value', strcmp(Plant.optimoptions.tspacing, 'logarithm'));
set(handles.manual, 'value', strcmp(Plant.optimoptions.tspacing, 'manual'));
set(handles.sequential, 'value', Plant.optimoptions.sequential);
set(handles.simultaneous, 'value', ~Plant.optimoptions.sequential);
set(handles.excessHeat, 'value', Plant.optimoptions.excessHeat);

if ismember('E',Plant.optimoptions.Outputs)
    set(handles.checkboxElectric,'value',1)
else set(handles.checkboxElectric,'value',0)
end
if ismember('H',Plant.optimoptions.Outputs)
    set(handles.checkboxHeat,'Value',1)
else set(handles.checkboxHeat,'Value',0)
end
if ismember('C',Plant.optimoptions.Outputs)
    set(handles.checkboxCooling,'Value',1)
else set(handles.checkboxCooling,'Value',0)
end
if ismember('S',Plant.optimoptions.Outputs)
    set(handles.checkboxSteam,'Value',1) 
else set(handles.checkboxSteam,'Value',0)
end
guidata(handles.figure1, handles);


% --- Outputs from this function are returned to the command line.
function varargout = OptimOptions_OutputFcn(hObject, eventdata, handles)

% --- Executes on button press in OK.
function OK_Callback(hObject, eventdata, handles)
global Plant
Plant.optimoptions.Interval = str2double(get(handles.Interval, 'String'));
Plant.optimoptions.Horizon = str2double(get(handles.Horizon, 'String'));
Plant.optimoptions.Resolution = str2double(get(handles.Resolution, 'String'));
Plant.optimoptions.Topt = str2double(get(handles.Topt, 'String'));
Plant.optimoptions.Tmpc = str2double(get(handles.Tmpc, 'String'));
Plant.optimoptions.nsSmooth = str2double(get(handles.nsSmooth, 'String'));
Plant.optimoptions.scaletime = str2double(get(handles.scaletime, 'String'));
Plant.optimoptions.thresholdSteps = str2double(get(handles.thresholdSteps, 'String'));
if get(handles.excessHeat, 'Value') ==1
    Plant.optimoptions.excessHeat = 1;%excessHeat is 1 if you can produce escess heat
else Plant.optimoptions.excessHeat = 0;%it is zero if all heat produced must be used
end
if get(handles.fastsimulation,'Value')==1
    Plant.optimoptions.fastsimulation = 1;
else Plant.optimoptions.fastsimulation = 0;
end
if get(handles.sequential,'Value')==1
    %% Only valid question if a) there are chillers and electric generators & b) there is cold thermal storage so that the chiller dispatch affects the electric dispatch
    %% Combined has the drawback of only single value for chiller efficiency (like hratio) rather than a quadratic fit
    Plant.optimoptions.sequential = 1;
else Plant.optimoptions.sequential = 0;
end
%% also build the list of outputs to consider (i.e. E, H, S, H2, C for electricity, heat steam hydrogen, cooling..)
Plant.optimoptions.Outputs = {};
if get(handles.checkboxElectric,'Value')==1
    Plant.optimoptions.Outputs(end+1) ={'E'}; 
end
if get(handles.checkboxHeat,'Value')==1
    Plant.optimoptions.Outputs(end+1) ={'H'}; 
end
if get(handles.checkboxCooling,'Value')==1
    Plant.optimoptions.Outputs(end+1) ={'C'}; 
end
if get(handles.checkboxSteam,'Value')==1
    Plant.optimoptions.Outputs(end+1) ={'S'}; 
end
close(gcf)

% --- Executes on button press in Cancel.
function Cancel_Callback(hObject, eventdata, handles)
close(gcf)

% --- Executes when selected object is changed in changingtimesteps.
function changingtimesteps_SelectionChangeFcn(hObject, eventdata, handles)
global Plant
switch get(eventdata.NewValue,'Tag')
    case 'constant'
        Plant.optimoptions.tspacing = 'constant';
    case 'manual'
        Plant.optimoptions.tspacing = 'manual';
        prompt = {'Specify a vector of times out of one horizon for each timestep: (ex: .05,.1,.25,.75,1 would result in one timestep at 5hr, 10hr, 25hr, 75hr, and 100hr for a 100 hour horizon)'};
        dlg_title = 'Manual Timesteps';
        num_lines = 1;
        def_ans = {'0.0035,0.0070,0.0105,0.0140,0.0175,0.0210,0.0417,0.0625,0.0833,0.1042,0.125,0.1667,0.2083,0.25,0.3333,0.4167,0.5,0.5833,0.6667,0.75,0.8333,0.9167,1'};
        a = inputdlg(prompt,dlg_title,num_lines,def_ans);
        Plant.optimoptions.manualT = str2double(strsplit(a{:}, ','));
    case 'linear'
        Plant.optimoptions.tspacing = 'constant';
    case 'logarithm'
        Plant.optimoptions.tspacing = 'logarithm';
end

function Interval_Callback(hObject, eventdata, handles)
function Horizon_Callback(hObject, eventdata, handles)
function Resolution_Callback(hObject, eventdata, handles)
function Topt_Callback(hObject, eventdata, handles)
function Tmpc_Callback(hObject, eventdata, handles)
function nsSmooth_Callback(hObject, eventdata, handles)
function scaletime_Callback(hObject, eventdata, handles)
function thresholdSteps_Callback(hObject, eventdata, handles)
function constant_Callback(hObject, eventdata, handles)
function logarithm_Callback(hObject, eventdata, handles)
function manual_Callback(hObject, eventdata, handles)
function slowsimulation_Callback(hObject, eventdata, handles)
function fastsimulation_Callback(hObject, eventdata, handles)
function simultaneous_Callback(hObject, eventdata, handles)
function sequential_Callback(hObject, eventdtaa, handles)


function Interval_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Horizon_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Resolution_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Topt_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Tmpc_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function nsSmooth_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function scaletime_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function checkboxElectric_Callback(hObject, eventdata, handles)
function checkboxHeat_Callback(hObject, eventdata, handles)
function checkboxCooling_Callback(hObject, eventdata, handles)
function checkboxSteam_Callback(hObject, eventdata, handles)

function thresholdSteps_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
