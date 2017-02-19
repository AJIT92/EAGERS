function varargout = ViewResultsGHG(varargin)
% VIEWRESULTSGHG MATLAB code for ViewResultsGHG.fig
%      VIEWRESULTSGHG, by itself, creates a new VIEWRESULTSGHG or raises the existing
%      singleton*.
%
%      H = VIEWRESULTSGHG returns the handle to a new VIEWRESULTSGHG or the handle to
%      the existing singleton*.
%
%      VIEWRESULTSGHG('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VIEWRESULTSGHG.M with the given input arguments.
%
%      VIEWRESULTSGHG('Property','Value',...) creates a new VIEWRESULTSGHG or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ViewResultsGHG_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ViewResultsGHG_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ViewResultsGHG

% Last Modified by GUIDE v2.5 10-Dec-2013 10:35:45

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ViewResultsGHG_OpeningFcn, ...
                   'gui_OutputFcn',  @ViewResultsGHG_OutputFcn, ...
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


% --- Executes just before ViewResultsGHG is made visible.
function ViewResultsGHG_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ViewResultsGHG (see VARARGIN)

% Choose default command line output for ViewResultsGHG
currentProject=max(varargin{1},1);
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);


global Model_dir
% UIWAIT makes ViewResultsGHG wait for user response (see UIRESUME)
% uiwait(handles.figure1);

projdir=fullfile(Model_dir, 'project'); %strrep(which('NREL_FCModel.m'),'\main\NREL_FCModel.m','\project');
files=dir(fullfile(projdir,'*.mat'));
list=strrep({files.name},'.mat','');

set(handles.popupmenuResultsFile,'string',list,'value',currentProject)

%plot options
plotSubregion={
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
set(handles.popupmenuSubregion,'string',plotSubregion,'value',1)
plotAxes(handles);

% --- Outputs from this function are returned to the command line.
function varargout = ViewResultsGHG_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in popupmenuResultsFile.
function popupmenuResultsFile_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuResultsFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Model_dir Project
str=get(hObject,'string');
val=get(hObject,'value');
load(fullfile(Model_dir, 'project',str{val})); %strrep(which('NREL_FCModel.m'),'\main\NREL_FCModel.m','\project')
%plotAxes(handles);


% --- Executes during object creation, after setting all properties.
function popupmenuResultsFile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuResultsFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pushbuttonUpdate.
function pushbuttonUpdate_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonUpdate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
plotAxes(handles)


% --- Executes on selection change in popupmenuSubregion.
function popupmenuSubregion_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuSubregion (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenuSubregion contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenuSubregion


% --- Executes during object creation, after setting all properties.
function popupmenuSubregion_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuSubregion (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function plotAxes(handles)
Ts = 0.25;
HeatToCO2 = 2.2822;%2.2822 kWhr/lbs CO2
BoilerEffect = 0.9;
HeatToNOx = [4.6844e-004 1.6718e-004];%lbs NOx/kWhr of fuel [large(>100 MMBtu/hr) small(<100 MMBtu/hr)]
HeatToSO2 = 0.000002;%lbs SO2/MWhr

global Project Model_dir
RESULTS = Project.Result;
eventdata = [];
popupmenuSubregion_Callback(handles.popupmenuSubregion,eventdata,handles)

maxDemandH = max(Project.Building.DemandH)*0.003412;%max heat demand in MMBtu/hr FIX THIS!!!
if maxDemandH < 100
    HeatToNOx = HeatToNOx(1);
else
    HeatToNOx = HeatToNOx(2);
end

regionVal = get(handles.popupmenuSubregion,'value');
regionList = {'AKGD';
                'AKMS';
                'ERCT';
                'FRCC';
                'HIMS';
                'HIOA';
                'MROE';
                'MROW';
                'NYLI';
                'NEWE';
                'NYCW';
                'NYUP';
                'RFCE';
                'RFCM';
                'RFCW';
                'SRMW';
                'SRMV';
                'SRSO';
                'SRTV';
                'SRVC';
                'SPNO';
                'SPSO';
                'CAMX';
                'NWPP';
                'RMPA';
                'AZNM';
            };

region = char(regionList(regionVal));
CO2CHP = RESULTS.eOut.CO2Total;%CO2 from generation system
NOxCHP = RESULTS.eOut.NOxTotal;%NOx from generation system
SO2CHP = RESULTS.eOut.SO2Total;%CO2 from generation system

emission_dir=fullfile(Model_dir,'Emissions','Emission Data','emissionRates.mat');
load(emission_dir);

gridRateCO2 = emissionRates(regionVal,4);
gridRateNOx = emissionRates(regionVal,1);
gridRateSO2 = emissionRates(regionVal,3);

CO2base = sum(Project.Building.DemandE)*Ts/1000*gridRateCO2;
NOxbase = sum(Project.Building.DemandE)*Ts/1000*gridRateNOx;
SO2base = sum(Project.Building.DemandE)*Ts/1000*gridRateSO2;

CO2offset = RESULTS.eOut.Total_Electricity_produced_kWh/1000*gridRateCO2;
NOxoffset = RESULTS.eOut.Total_Electricity_produced_kWh/1000*gridRateNOx;
SO2offset = RESULTS.eOut.Total_Electricity_produced_kWh/1000*gridRateSO2;

TotalHeatkWh = RESULTS.eOut.Total_DG_Heat_out_kWh + RESULTS.eOut.Total_peak_burner_heat_out_kWh;
BaseBoilerCO2 = TotalHeatkWh/BoilerEffect/HeatToCO2;
NewBoilerCO2 =  RESULTS.eOut.Total_peak_burner_heat_out_kWh/BoilerEffect/HeatToCO2;
BaseBoilerNOx = TotalHeatkWh/BoilerEffect*HeatToNOx;
NewBoilerNOx =  RESULTS.eOut.Total_peak_burner_heat_out_kWh/BoilerEffect*HeatToNOx;
BaseBoilerSO2 = TotalHeatkWh/BoilerEffect*HeatToSO2;
NewBoilerSO2 =  RESULTS.eOut.Total_peak_burner_heat_out_kWh/BoilerEffect*HeatToSO2;

%%Plot CO2
CO2new = [(CO2base-CO2offset) CO2CHP NewBoilerCO2];
Fig.Y1 = [CO2base 0 BaseBoilerCO2; CO2new];
axes(handles.axes1)
bar(Fig.Y1,'stacked')
set(gca,'XTickLabel',{'Baseline Emissions', 'CHP Emissions'})
ylabel('lbs CO2')
l = legend('Grid emissions', 'CHP emissions','Boiler emissions');
set(l, 'fontsize', 8);

%%Plot NOx
NOxnew = [(NOxbase-NOxoffset) NOxCHP NewBoilerNOx];
Fig.Y2 = [NOxbase 0 BaseBoilerNOx; NOxnew];
axes(handles.axes2)
bar(Fig.Y2,'stacked')
set(gca,'XTickLabel',{'Baseline Emissions', 'CHP Emissions'})
ylabel('lbs NOx')
l = legend('Grid emissions', 'CHP emissions', 'Boiler emissions');
set(l, 'fontsize', 8);
%%Plot SO2
SO2new = [(SO2base-SO2offset) SO2CHP NewBoilerSO2];
Fig.Y3 = [SO2base 0 BaseBoilerSO2; SO2new];
axes(handles.axes3)
bar(Fig.Y3,'stacked')
set(gca,'XTickLabel',{'Baseline Emissions', 'CHP Emissions'})
ylabel('lbs SO2')
l = legend('Grid emissions', 'CHP emissions','Boiler emissions');
set(l, 'fontsize', 8);