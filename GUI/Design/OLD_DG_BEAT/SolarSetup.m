function varargout = SolarSetup(varargin)
% SOLARSETUP MATLAB code for SolarSetup.fig
%      SOLARSETUP, by itself, creates a new SOLARSETUP or raises the existing
%      singleton*.
%
%      H = SOLARSETUP returns the handle to a new SOLARSETUP or the handle to
%      the existing singleton*.
%
%      SOLARSETUP('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SOLARSETUP.M with the given input arguments.
%
%      SOLARSETUP('Property','Value',...) creates a new SOLARSETUP or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SolarSetup_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SolarSetup_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SolarSetup

% Last Modified by GUIDE v2.5 18-Apr-2013 14:46:16

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SolarSetup_OpeningFcn, ...
                   'gui_OutputFcn',  @SolarSetup_OutputFcn, ...
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


% --- Executes just before SolarSetup is made visible.
function SolarSetup_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Add a wait dialogue box
figure_Sol = get(0,'CurrentFigure');
MSG_Sol=msgbox('Loading');
figure(figure_Sol)

% varargin   command line arguments to SolarSetup (see VARARGIN)
global Project SYSINDEX  SolarOriginal
SolarOriginal = Project.Renewable.Solar(SYSINDEX);
% Choose default command line output for SolarSetup
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
solar = Project.Renewable.Solar(SYSINDEX);
stateName = {'Alabama';'Alaska';'Arizona';'Arkansas';'California';'Colorado';'Connecticut';'Delaware';'Florida';'Georgia';'Hawaii';'Idaho';'Illinois';'Indiana';'Iowa';'Kansas';
             'Kentucky';'Louisiana';'Maine';'Maryland';'Massachusetts';'Michigan';'Minnesota';'Mississippi';'Missouri';'Montana';'Nebraska';'Nevada';'New Hampshire';'New Jersey';
             'New Mexico';'New York';'North Carolina';'North Dakota';'Ohio';'Oklahoma';'Oregon';'Pennsylvania';'Rhode Island';'South Carolina';'South Dakota';'Tennessee';'Texas';
             'Utah';'Vermont';'Virginia';'Washington';'West Virginia';'Wisconsin';'Wyoming';};
stateNum = find(strcmp(solar.State,stateName));

set(handles.popupmenuState,'string',stateName,'value',stateNum)        
set(handles.editName,'string',solar.Name)
set(handles.uitableDCAC,'Data',solar.Data)
set(handles.editSize,'string',solar.Size)
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

solar = Project.Renewable.Solar;
CHPsize = 0;
for i = 1:1:length(Project.System.CHP)
    CHPsize = CHPsize + Project.System.CHP(i).SysSize(1);
end
AnnualGen = 8760*CHPsize;
AnnualDem = sum(Project.Building.DemandE)*8760/length(Project.Building.DemandE);
AnnualPower = zeros(length(solar),1);
for i = 1:1:length(solar)
    AnnualPower(i) = annualPower(Project.Renewable.Solar(i));
    solar(i).GenFrac = AnnualPower(i)/AnnualGen*100;
    solar(i).DemFrac = AnnualPower(i)/AnnualDem*100;
end
set(handles.editSizeAgen,'string',solar(SYSINDEX).GenFrac)
set(handles.editSizeAdem,'string',solar(SYSINDEX).DemFrac)
Project.Renewable.Solar = solar;

close(MSG_Sol)

% --- Outputs from this function are returned to the command line.
function varargout = SolarSetup_OutputFcn(hObject, eventdata, handles) 
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
global Project SYSINDEX
Project.Renewable.Solar(SYSINDEX).Name=get(hObject,'string');

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
global Project SYSINDEX Model_dir
val=get(hObject,'value');
str=get(hObject,'string');
Project.Renewable.Solar(SYSINDEX).State=str{val};
SolarDir = fullfile(Model_dir, 'System Library','Solar'); %strrep(which('NREL_FCModel.m'),'\main\NREL_FCModel.m','\component library\Solar');
if strcmp(Project.Renewable.Solar(SYSINDEX).PVtype,'flat')
    load(strcat(SolarDir,[filesep 'solarData' filesep 'GlobalHorizontal.mat']));
    Project.Renewable.Solar(SYSINDEX).Irrad = GlobalHorizontal(:,val);
else load(strcat(SolarDir,[filesep 'solarData' filesep 'DirectNormal.mat']));
    Project.Renewable.Solar(SYSINDEX).Irrad = DirectNormal(:,val);
end
load(strcat(SolarDir,[filesep 'solarData' filesep 'Azimuth.mat']));
load(strcat(SolarDir,[filesep 'solarData' filesep 'Zenith.mat']));
Project.Renewable.Solar(SYSINDEX).SunAz = Azimuth(:,val);
Project.Renewable.Solar(SYSINDEX).SunZen = Zenith(:,val);

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
global Project SYSINDEX
oldGenFrac = str2double(get(handles.editSizeAgen,'string'))/100;
oldDemFrac = str2double(get(handles.editSizeAdem,'string'))/100;
oldSize = Project.Renewable.Solar(SYSINDEX).Size;
oldSizem2 = Project.Renewable.Solar(SYSINDEX).Sizem2;
newSize = str2double(get(hObject,'string'));
Project.Renewable.Solar(SYSINDEX).Size= newSize;
Project.Renewable.Solar(SYSINDEX).Sizem2= oldSizem2*(newSize/oldSize);
set(handles.editSize,'string',Project.Renewable.Solar(SYSINDEX).Size)
set(handles.editSizem2,'string',Project.Renewable.Solar(SYSINDEX).Sizem2)
set(handles.editSizeAgen,'string',oldGenFrac*(newSize/oldSize)*100)
set(handles.editSizeAdem,'string',oldDemFrac*(newSize/oldSize)*100)

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
global Project SYSINDEX
Project.Renewable.Solar(SYSINDEX).Eff=str2double(get(hObject,'string'));
set(handles.editEfficiency,'string',Project.Renewable.Solar(SYSINDEX).Eff)
oldGenFrac = str2double(get(handles.editSizeAgen,'string'))/100;
oldDemFrac = str2double(get(handles.editSizeAdem,'string'))/100;
oldSize = Project.Renewable.Solar(SYSINDEX).Size;
newSize = Project.Renewable.Solar(SYSINDEX).Sizem2*Project.Renewable.Solar(SYSINDEX).Eff;
Project.Renewable.Solar(SYSINDEX).Size= newSize;
set(handles.editSize,'string',Project.Renewable.Solar(SYSINDEX).Size)
set(handles.editSizeAgen,'string',oldGenFrac*(newSize/oldSize)*100)
set(handles.editSizeAdem,'string',oldDemFrac*(newSize/oldSize)*100)

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
global Project SYSINDEX
oldGenFrac = str2double(get(handles.editSizeAgen,'string'))/100;
oldDemFrac = str2double(get(handles.editSizeAdem,'string'))/100;
oldSize = Project.Renewable.Solar(SYSINDEX).Size;
oldSizem2 = Project.Renewable.Solar(SYSINDEX).Sizem2;
newSizem2 = str2double(get(hObject,'string'));
Project.Renewable.Solar(SYSINDEX).Sizem2= newSize;
Project.Renewable.Solar(SYSINDEX).Size= oldSize*(newSizem2/oldSizem2);
set(handles.editSize,'string',Project.Renewable.Solar(SYSINDEX).Size)
set(handles.editSizem2,'string',Project.Renewable.Solar(SYSINDEX).Sizem2)
set(handles.editSizeAgen,'string',oldGenFrac*(newSizem2/oldSizem2)*100)
set(handles.editSizeAdem,'string',oldDemFrac*(newSizem2/oldSizem2)*100)

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
global Project SYSINDEX
Project.Renewable.Solar(SYSINDEX).Tilt=str2double(get(hObject,'string'));
set(handles.editTilt,'string',Project.Renewable.Solar(SYSINDEX).Tilt)

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
global Project SYSINDEX
Project.Renewable.Solar(SYSINDEX).Azimuth=str2double(get(hObject,'string'));
set(handles.editAzimuth,'string',Project.Renewable.Solar(SYSINDEX).Azimuth)

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
global Project SYSINDEX
if strcmp(Project.Renewable.Solar(SYSINDEX).Tracking,'fixed')
    track = 1;
elseif strcmp(Project.Renewable.Solar(SYSINDEX).Tracking,'1axis')
    track = 2;
else track =3;
end
Project.Renewable.Solar(SYSINDEX).Tilt = OptimalZenith(Project.Renewable.Solar(SYSINDEX).Irrad,Project.Renewable.Solar(SYSINDEX).SunAz,Project.Renewable.Solar(SYSINDEX).SunZen,track,Project.Renewable.Solar(SYSINDEX).Azimuth,Project.Renewable.Solar(SYSINDEX).Eff);
set(handles.editTilt,'string',Project.Renewable.Solar(SYSINDEX).Tilt)

% --- Executes on button press in pushbuttonSaveOnlyToProject.
function pushbuttonSaveOnlyToProject_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonSaveOnlyToProject (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clear global SYSINDEX SolarOriginal
close(gcf)
uiresume

% --- Executes on button press in pushbuttonCancel.
function pushbuttonCancel_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonCancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Project SYSINDEX SolarOriginal
Project.Renewable.Solar(SYSINDEX) = SolarOriginal;
clear global SYSINDEX SolarOriginal
close(gcf)


% --- Executes on button press in pushbuttonSaveAs.
function pushbuttonSaveAs_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonSaveAs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Project SYSINDEX Model_dir
[f,p]=uiputfile('*.mat','Save As Solar Component',fullfile(Model_dir,'System Library','Solar'));
if f==0;return;end
component=Project.Renewable.Solar(SYSINDEX);
save([p f],'component')
clear global SYSINDEX SolarOriginal 
close(gcf)
uiresume


% --- Executes when selected object is changed in uipanelType.
function uipanelType_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanelType 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
global Project SYSINDEX
switch get(eventdata.NewValue,'Tag')
    case 'radiobuttonFixed'
        Project.Renewable.Solar(SYSINDEX).PVtype = 'flat';
    case 'radiobuttonConcentrated'
        Project.Renewable.Solar(SYSINDEX).PVtype = 'concentrated';
end


% --- Executes when selected object is changed in uipanelTrack.
function uipanelTrack_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanelTrack 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
global Project SYSINDEX
switch get(eventdata.NewValue,'Tag')
    case 'radiobuttonFixed'
        Project.Renewable.Solar(SYSINDEX).Tracking = 'fixed';
    case 'radiobutton1axis'
        Project.Renewable.Solar(SYSINDEX).Tracking = '1axis';
    case 'radiobutton2axis'
        Project.Renewable.Solar(SYSINDEX).Tracking = '2axis';
end



function editSizeAgen_Callback(hObject, eventdata, handles)
% hObject    handle to editSizeAgen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Project SYSINDEX
SizeAgen = str2double(get(hObject,'String'));
SolarSizing(SYSINDEX,'generation','no',SizeAgen)
set(handles.editSize,'string',Project.Renewable.Solar(SYSINDEX).Size)
set(handles.editSizem2,'string',Project.Renewable.Solar(SYSINDEX).Sizem2)
set(handles.editSizeAgen,'string',Project.Renewable.Solar(SYSINDEX).GenFrac)
set(handles.editSizeAdem,'string',Project.Renewable.Solar(SYSINDEX).DemFrac)

% --- Executes during object creation, after setting all properties.
function editSizeAgen_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editSizeAgen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editSizeAdem_Callback(hObject, eventdata, handles)
% hObject    handle to editSizeAdem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Project SYSINDEX
SizeAdem = str2double(get(hObject,'String'));
SolarSizing(SYSINDEX,'demand','no',SizeAdem)
set(handles.editSize,'string',Project.Renewable.Solar(SYSINDEX).Size)
set(handles.editSizem2,'string',Project.Renewable.Solar(SYSINDEX).Sizem2)
set(handles.editSizeAgen,'string',Project.Renewable.Solar(SYSINDEX).GenFrac)
set(handles.editSizeAdem,'string',Project.Renewable.Solar(SYSINDEX).DemFrac)

% --- Executes during object creation, after setting all properties.
function editSizeAdem_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editSizeAdem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
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
