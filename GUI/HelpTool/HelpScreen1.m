function varargout = HelpScreen1(varargin)
% HELPSCREEN1 MATLAB code for HelpScreen1.fig
%      HELPSCREEN1, by itself, creates a new HELPSCREEN1 or raises the existing
%      singleton*.
%
%      H = HELPSCREEN1 returns the handle to a new HELPSCREEN1 or the handle to
%      the existing singleton*.
%
%      HELPSCREEN1('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in HELPSCREEN1.M with the given input arguments.
%
%      HELPSCREEN1('Property','Value',...) creates a new HELPSCREEN1 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before HelpScreen1_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to HelpScreen1_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help HelpScreen1

% Last Modified by GUIDE v2.5 27-Jan-2014 17:02:07

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @HelpScreen1_OpeningFcn, ...
                   'gui_OutputFcn',  @HelpScreen1_OutputFcn, ...
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


% --- Executes just before HelpScreen1 is made visible.
function HelpScreen1_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to HelpScreen1 (see VARARGIN)

% Choose default command line output for HelpScreen1
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
% global Model_dir
% Pic_dir=fullfile(Model_dir, 'graphics');
% pFileName = fullfile(Pic_dir,'nrel_logo_large.jpg');
% h=handles.axes1;
% image(imread(pFileName),'Parent',h)
% axis(h,'image')
% axis(h,'off')


% --- Outputs from this function are returned to the command line.
function varargout = HelpScreen1_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbuttonExplore.
function pushbuttonExplore_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonExplore (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close(gcf)
HelpScreen2a()

% --- Executes on button press in pushbuttonSingleBuild.
function pushbuttonSingleBuild_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonSingleBuild (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close(gcf)
HelpScreen2b()

% --- Executes on button press in pushbuttonSingleFC.
function pushbuttonSingleFC_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonSingleFC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close(gcf)
HelpScreen2c()

% --- Executes on button press in pushbuttonFleet.
function pushbuttonFleet_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonFleet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close(gcf)
HelpScreen2d()

% --- Executes on button press in pushbuttonNatSurvey.
function pushbuttonNatSurvey_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonNatSurvey (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close(gcf)
HelpScreen2e()
