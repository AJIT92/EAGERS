function varargout = HelpScreen4b(varargin)
% HELPSCREEN4B MATLAB code for HelpScreen4b.fig
%      HELPSCREEN4B, by itself, creates a new HELPSCREEN4B or raises the existing
%      singleton*.
%
%      H = HELPSCREEN4B returns the handle to a new HELPSCREEN4B or the handle to
%      the existing singleton*.
%
%      HELPSCREEN4B('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in HELPSCREEN4B.M with the given input arguments.
%
%      HELPSCREEN4B('Property','Value',...) creates a new HELPSCREEN4B or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before HelpScreen4b_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to HelpScreen4b_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help HelpScreen4b

% Last Modified by GUIDE v2.5 27-Jan-2014 12:39:45

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @HelpScreen4b_OpeningFcn, ...
                   'gui_OutputFcn',  @HelpScreen4b_OutputFcn, ...
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


% --- Executes just before HelpScreen4b is made visible.
function HelpScreen4b_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to HelpScreen4b (see VARARGIN)

% Choose default command line output for HelpScreen4b
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
global Model_dir
Pic_dir=fullfile(Model_dir, 'GUI','HelpTool','screenshots');
pFileName = fullfile(Pic_dir,'HelpScreen4b.jpg');
h=handles.axes1;
image(imread(pFileName),'Parent',h)
axis(h,'image')
axis(h,'off');
pFileName = fullfile(Pic_dir,'HelpScreen4b2.jpg');
h=handles.axes2;
image(imread(pFileName),'Parent',h)
axis(h,'image')
axis(h,'off');


% --- Outputs from this function are returned to the command line.
function varargout = HelpScreen4b_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbuttonNext.
function pushbuttonNext_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonNext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close(gcf)
HelpScreen5b()
