function varargout = SetupSolar(varargin)
% SETUPSOLAR MATLAB code for SetupSolar.fig
%      SETUPSOLAR, by itself, creates a new SETUPSOLAR or raises the existing
%      singleton*.
%
%      H = SETUPSOLAR returns the handle to a new SETUPSOLAR or the handle to
%      the existing singleton*.
%
%      SETUPSOLAR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SETUPSOLAR.M with the given input arguments.
%
%      SETUPSOLAR('Property','Value',...) creates a new SETUPSOLAR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SetupSolar_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SetupSolar_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SetupSolar

% Last Modified by GUIDE v2.5 26-Nov-2014 13:34:00

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SetupSolar_OpeningFcn, ...
                   'gui_OutputFcn',  @SetupSolar_OutputFcn, ...
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


% --- Executes just before SetupSolar is made visible.
function SetupSolar_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SetupSolar (see VARARGIN)
global Plant SYSINDEX  Original figHandle
figHandle = handles.figure1;
Original = Plant.Generator(SYSINDEX);
% Choose default command line output for SetupSolar
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
solar = Plant.Generator(SYSINDEX).VariableStruct;
stateName = {'Alabama';'Alaska';'Arizona';'Arkansas';'California';'Colorado';'Connecticut';'Delaware';'Florida';'Georgia';'Hawaii';'Idaho';'Illinois';'Indiana';'Iowa';'Kansas';
             'Kentucky';'Louisiana';'Maine';'Maryland';'Massachusetts';'Michigan';'Minnesota';'Mississippi';'Missouri';'Montana';'Nebraska';'Nevada';'New Hampshire';'New Jersey';
             'New Mexico';'New York';'North Carolina';'North Dakota';'Ohio';'Oklahoma';'Oregon';'Pennsylvania';'Rhode Island';'South Carolina';'South Dakota';'Tennessee';'Texas';
             'Utah';'Vermont';'Virginia';'Washington';'West Virginia';'Wisconsin';'Wyoming';};
stateNum = find(strcmp(solar.State,stateName));

set(handles.popupmenuState,'string',stateName,'value',stateNum)        
set(handles.editName,'string',Plant.Generator(SYSINDEX).Name)
set(handles.uitableDCAC,'Data',solar.Data)
set(handles.editSize,'string',Plant.Generator(SYSINDEX).Size)
set(handles.editEfficiency,'string',solar.Eff)
set(handles.editSizem2,'string',solar.Sizem2)
set(handles.editTilt,'string',solar.Tilt)
set(handles.editAzimuth,'string',solar.Azimuth)
if strcmp(solar.PVtype, 'flat')
    set(handles.uipanelType,'SelectedObject', handles.radiobuttonFlat)
else    set(handles.uipanelType,'SelectedObject', handles.radiobuttonConcentrated)
end
if strcmp(solar.Tracking, 'fixed')
    set(handles.uipanelTrack,'SelectedObject', handles.radiobuttonFixed)
elseif strcmp(solar.Tracking, '1axis')
    set(handles.uipanelTrack,'SelectedObject', handles.radiobutton1axis)
else set(handles.uipanelTrack,'SelectedObject', handles.radiobutton2axis)
end

nSys = length(Plant.Generator);
nSolar = 0;
NGsize = 0;
for g = 1:1:nSys
    if strcmp(Plant.Generator(g).Type,'Renewable') && strcmp(Plant.Generator(g).Source,'Solar')
        nSolar = nSolar+1;
        AnnualPower(nSolar) = annualPower(Plant.Generator(g).VariableStruct);
    elseif strcmp(Plant.Generator(g).Source,'NG')
        NGsize = NGsize + Plant.Generator(g).Size;
    end
end
nSolar = 0;
AnnualGen = 8760*NGsize;
AnnualDem = mean(Plant.Data.Demand.E*1000)*365*24;
for g = 1:1:nSys
    if strcmp(Plant.Generator(g).Type,'Renewable') && strcmp(Plant.Generator(g).Source,'Solar')
        nSolar = nSolar+1;
        Plant.Generator(g).VariableStruct.GenFrac = AnnualPower(nSolar)/AnnualGen*100;
        Plant.Generator(g).VariableStruct.DemFrac = AnnualPower(nSolar)/AnnualDem*100;
    end
end


% --- Outputs from this function are returned to the command line.
function varargout = SetupSolar_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function editName_Callback(hObject, eventdata, handles)
% hObject    handle to editName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Plant SYSINDEX
Plant.Generator(SYSINDEX).Name=get(hObject,'string');

% --- Executes during object creation, after setting all properties.
function editName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in popupmenuState.
function popupmenuState_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuState (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Plant SYSINDEX Model_dir
val=get(hObject,'value');
str=get(hObject,'string');
Plant.Generator(SYSINDEX).VariableStruct.State=str{val};
SolarDir = fullfile(Model_dir, '\component library\Solar'); %strrep(which('NREL_FCModel.m'),'\main\NREL_FCModel.m','\component library\Solar');
if strcmp(Plant.Generator(SYSINDEX).VariableStruct.PVtype,'flat')
    load(strcat(SolarDir,'\solarData\GlobalHorizontal.mat'));
    Plant.Generator(SYSINDEX).VariableStruct.Irrad = GlobalHorizontal(:,val);
else load(strcat(SolarDir,'\solarData\DirectNormal.mat'));
    Plant.Generator(SYSINDEX).VariableStruct.Irrad = DirectNormal(:,val);
end
load(strcat(SolarDir,'\solarData\Azimuth.mat'));
load(strcat(SolarDir,'\solarData\Zenith.mat'));
Plant.Generator(SYSINDEX).VariableStruct.SunAz = Azimuth(:,val);
Plant.Generator(SYSINDEX).VariableStruct.SunZen = Zenith(:,val);

% --- Executes during object creation, after setting all properties.
function popupmenuState_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuState (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editSize_Callback(hObject, eventdata, handles)
% hObject    handle to editSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Plant SYSINDEX
Plant.Generator(SYSINDEX).Size= str2double(get(hObject,'string'));
Plant.Generator(SYSINDEX).VariableStruct.Sizem2= Plant.Generator(SYSINDEX).Size/Plant.Generator(SYSINDEX).VariableStruct.Eff;
set(handles.editSizem2,'string',Plant.Generator(SYSINDEX).VariableStruct.Sizem2)

% --- Executes during object creation, after setting all properties.
function editSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editEfficiency_Callback(hObject, eventdata, handles)
% hObject    handle to editEfficiency (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Plant SYSINDEX
Plant.Generator(SYSINDEX).VariableStruct.Eff=str2double(get(hObject,'string'));
Plant.Generator(SYSINDEX).VariableStruct.Sizem2= Plant.Generator(SYSINDEX).Size/Plant.Generator(SYSINDEX).VariableStruct.Eff;
set(handles.editSizem2,'string',Plant.Generator(SYSINDEX).VariableStruct.Sizem2)

% --- Executes during object creation, after setting all properties.
function editEfficiency_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editEfficiency (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editSizem2_Callback(hObject, eventdata, handles)
% hObject    handle to editSizem2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Plant SYSINDEX
Plant.Generator(SYSINDEX).VariableStruct.Sizem2= str2double(get(hObject,'string'));
Plant.Generator(SYSINDEX).Size= Plant.Generator(SYSINDEX).VariableStruct.Sizem2*Plant.Generator(SYSINDEX).VariableStruct.Eff;
set(handles.editSize,'string',Plant.Generator(SYSINDEX).Size)

% --- Executes during object creation, after setting all properties.
function editSizem2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editSizem2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editTilt_Callback(hObject, eventdata, handles)
% hObject    handle to editTilt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Plant SYSINDEX
Plant.Generator(SYSINDEX).VariableStruct.Tilt=str2double(get(hObject,'string'));
set(handles.editTilt,'string',Plant.Generator(SYSINDEX).VariableStruct.Tilt)

% --- Executes during object creation, after setting all properties.
function editTilt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editTilt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editAzimuth_Callback(hObject, eventdata, handles)
% hObject    handle to editAzimuth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Plant SYSINDEX
Plant.Generator(SYSINDEX).VariableStruct.Azimuth=str2double(get(hObject,'string'));
set(handles.editAzimuth,'string',Plant.Generator(SYSINDEX).VariableStruct.Azimuth)

% --- Executes during object creation, after setting all properties.
function editAzimuth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editAzimuth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbuttonOptZenith.
function pushbuttonOptZenith_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonOptZenith (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Plant SYSINDEX
if strcmp(Plant.Generator(SYSINDEX).VariableStruct.Tracking,'fixed')
    track = 1;
elseif strcmp(Plant.Generator(SYSINDEX).VariableStruct.Tracking,'1axis')
    track = 2;
else track =3;
end
Plant.Generator(SYSINDEX).VariableStruct.Tilt = OptimalZenith(Plant.Generator(SYSINDEX).VariableStruct.Irrad,Plant.Generator(SYSINDEX).VariableStruct.SunAz,Plant.Generator(SYSINDEX).VariableStruct.SunZen,track,Plant.Generator(SYSINDEX).VariableStruct.Azimuth,Plant.Generator(SYSINDEX).VariableStruct.Eff);
set(handles.editTilt,'string',Plant.Generator(SYSINDEX).VariableStruct.Tilt)

% --- Executes on button press in pushbuttonSaveOnlyToPlant.
function pushbuttonSaveOnlyToPlant_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonSaveOnlyToPlant (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close(gcf)

% --- Executes on button press in pushbuttonCancel.
function pushbuttonCancel_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonCancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Plant SYSINDEX Original
Plant.Generator(SYSINDEX) = Original;
clear global Original
close(gcf)

% --- Executes on button press in pushbuttonSaveAs.
function pushbuttonSaveAs_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonSaveAs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Plant SYSINDEX Model_dir
savedir= fullfile(Model_dir, 'component library','Solar',strcat(Plant.Generator(SYSINDEX).Name,'.mat')); 
[f,p]=uiputfile(savedir,'Save As Solar Component');
if f==0;return;end
if isfield(Plant.Generator(SYSINDEX),'OpMatA')%if the model has already been run
    component = rmfield(Plant.Generator(SYSINDEX),'OpMatA');%remove OpMatA when saving the new component
else
    component=Plant.Generator(SYSINDEX);    
end
save([p f],'component')
clear global Original
close(gcf)


% --- Executes when selected object is changed in uipanelType.
function uipanelType_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanelType 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
global Plant SYSINDEX
switch get(eventdata.NewValue,'Tag')
    case 'radiobuttonFixed'
        Plant.Generator(SYSINDEX).VariableStruct.PVtype = 'flat';
    case 'radiobuttonConcentrated'
        Plant.Generator(SYSINDEX).VariableStruct.PVtype = 'concentrated';
end


% --- Executes when selected object is changed in uipanelTrack.
function uipanelTrack_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanelTrack 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
global Plant SYSINDEX
switch get(eventdata.NewValue,'Tag')
    case 'radiobuttonFixed'
        Plant.Generator(SYSINDEX).VariableStruct.Tracking = 'fixed';
    case 'radiobutton1axis'
        Plant.Generator(SYSINDEX).VariableStruct.Tracking = '1axis';
    case 'radiobutton2axis'
        Plant.Generator(SYSINDEX).VariableStruct.Tracking = '2axis';
end


function Tilt = OptimalZenith(Irrad,Az,Zen,track,Ao,Eff)
minZen = min(Zen+99*(Zen==0));
meanZen = sum(Zen)/sum(Zen>0);

% Test multiple tilt angles
To = linspace(minZen,meanZen,20);
Col = linspace(1,1,length(To));
if track ==3
    Tilt = 0;
elseif track ==1
    Po = (Irrad*Col).*cosd(Zen*Col-linspace(1,1,length(Zen))'*To).*(cosd(Az-Ao)*Col)*Eff;
    AnPow = sum(Po,1);
    [MaxPow, I] = max(AnPow);
    if I ==1
        I =2;
    end
    To = linspace(To(I-1),To(I+1),20);
    Po = (Irrad*Col).*cosd(Zen*Col-linspace(1,1,length(Zen))'*To).*(cosd(Az-Ao)*Col)*Eff;
    AnPow = sum(Po,1);
    [MaxPow, I] = max(AnPow);
    Tilt = To(I);
elseif track ==2
    Po = (Irrad*Col).*cosd(Zen*Col-linspace(1,1,length(Zen))'*To)*Eff;
    AnPow = sum(Po,1);
    [MaxPow, I] = max(AnPow);
    if I ==1
        I =2;
    end
    To = linspace(To(I-1),To(I+1),20);
    Po = (Irrad*Col).*cosd(Zen*Col-linspace(1,1,length(Zen))'*To).*(cosd(Az-Ao)*Col)*Eff;
    AnPow = sum(Po,1);
    [MaxPow, I] = max(AnPow);
    Tilt = To(I);
end

function Power = annualPower(solar)
if strcmp(solar.Tracking,'fixed')
    Power = 8760*solar.Sizem2*sum((solar.Irrad/1000).*cosd(solar.SunZen-solar.Tilt).*(cosd(solar.SunAz-solar.Azimuth))*solar.Eff)/length(solar.SunAz);
elseif strcmp(solar.Tracking,'1axis')
    Power = 8760*solar.Sizem2*sum((solar.Irrad/1000).*cosd(solar.SunZen-solar.Tilt)*solar.Eff)/length(solar.SunAz);
else Power = 8760*solar.Sizem2*sum((solar.Irrad/1000)*solar.Eff)/length(solar.SunAz);
end
