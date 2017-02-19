function varargout = ThermalStorageSetup(varargin)
% THERMALSTORAGESETUP M-file for ThermalStorageSetup.fig
%      THERMALSTORAGESETUP, by itself, creates a new THERMALSTORAGESETUP or raises the existing
%      singleton*.
%
%      H = THERMALSTORAGESETUP returns the handle to a new THERMALSTORAGESETUP or the handle to
%      the existing singleton*.
%
%      THERMALSTORAGESETUP('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in THERMALSTORAGESETUP.M with the given input arguments.
%
%      THERMALSTORAGESETUP('Property','Value',...) creates a new THERMALSTORAGESETUP or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ThermalStorageSetup_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ThermalStorageSetup_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ThermalStorageSetup

% Last Modified by GUIDE v2.5 21-Feb-2013 09:44:06

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ThermalStorageSetup_OpeningFcn, ...
                   'gui_OutputFcn',  @ThermalStorageSetup_OutputFcn, ...
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


% --- Executes just before ThermalStorageSetup is made visible.
function ThermalStorageSetup_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ThermalStorageSetup (see VARARGIN)

% Add a wait dialogue box
figure_TS = get(0,'CurrentFigure');
MSG_TS=msgbox('Loading');
figure(figure_TS)

% Choose default command line output for ThermalStorageSetup
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

global Project SYSINDEX TESoriginal
TESoriginal = Project.System.TES;
OptimalTESsizing('editSize');
COMPONENT=Project.System.TES(SYSINDEX);

set(handles.editName,'string',COMPONENT.Name)
set(handles.textType,'string',COMPONENT.Type)
set(handles.editVersion,'string',COMPONENT.Version)
set(handles.editDescription,'string',COMPONENT.Description)
% set(handles.editPicture,'string',COMPONENT.Picture)
% editPicture_Callback(handles.editPicture,eventdata,handles)
set(handles.editTcold,'string',COMPONENT.Tcold)
set(handles.editThot,'string',COMPONENT.Thot)
set(handles.editLossDay,'string',COMPONENT.LossDay)
set(handles.editFillRate,'string',COMPONENT.FillRate)
set(handles.editFillRatePerc,'string',COMPONENT.FillRatePerc)
set(handles.editDischRate,'string',COMPONENT.DischRate)
set(handles.editDischRatePerc,'string',COMPONENT.DischRatePerc)
ButtonMenu = handles.uipanelTEStype;
radio1 = handles.radiobuttonColdTES;
radio2 = handles.radiobuttonHotTES;
if strcmp(COMPONENT.TEStype, 'cold')
    set(ButtonMenu,'SelectedObject', radio1)
else
    set(ButtonMenu,'SelectedObject', radio2)
end
set(handles.editSizeGal,'string',COMPONENT.SizeGal)
set(handles.editSize,'string',COMPONENT.SysSize)
set(handles.editOptSize,'string',COMPONENT.OptSize)
set(handles.editSizeHours,'string',COMPONENT.SizeHours)
fillRatesStatic(COMPONENT,handles)

close(MSG_TS)

% --- Outputs from this function are returned to the command line.
function varargout = ThermalStorageSetup_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function editSize_Callback(hObject, eventdata, handles)
% hObject    handle to editSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editSize as text
%        str2double(get(hObject,'String')) returns contents of editSize as a double
global Project SYSINDEX
Project.System.TES(SYSINDEX).SysSize=str2double(get(hObject,'string'));
OptimalTESsizing('editSize');
COMPONENT = Project.System.TES(SYSINDEX);
set(handles.editSizeGal,'string',COMPONENT.SizeGal)
set(handles.editSize,'string',COMPONENT.SysSize)
set(handles.editOptSize,'string',COMPONENT.OptSize)
set(handles.editSizeHours,'string',COMPONENT.SizeHours)
fillRatesStatic(COMPONENT,handles)


% --- Executes during object creation, after setting all properties.
function editSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function editOptSize_Callback(hObject, eventdata, handles)
% hObject    handle to editOptSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Project SYSINDEX
Project.System.TES(SYSINDEX).OptSize=str2double(get(hObject,'string'));
OptimalTESsizing('OptSize');
COMPONENT = Project.System.TES(SYSINDEX);
set(handles.editSizeGal,'string',COMPONENT.SizeGal)
set(handles.editSize,'string',COMPONENT.SysSize)
set(handles.editOptSize,'string',COMPONENT.OptSize)
set(handles.editSizeHours,'string',COMPONENT.SizeHours)
fillRatesStatic(COMPONENT,handles)

% --- Executes during object creation, after setting all properties.
function editOptSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editOptSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function editName_Callback(hObject, eventdata, handles)
% hObject    handle to editName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Project SYSINDEX
Project.System.TES(SYSINDEX).Name=get(hObject,'string');

% --- Executes during object creation, after setting all properties.
function editName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editDescription_Callback(hObject, eventdata, handles)
% hObject    handle to editDescription (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Project SYSINDEX
Project.System.TES(SYSINDEX).Description=get(hObject,'string');

% --- Executes during object creation, after setting all properties.
function editDescription_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editDescription (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editVersion_Callback(hObject, eventdata, handles)
% hObject    handle to editVersion (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Project SYSINDEX
Project.System.TES(SYSINDEX).Version=get(hObject,'string');

% --- Executes during object creation, after setting all properties.
function editVersion_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editVersion (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editPicture_Callback(hObject, eventdata, handles)
% hObject    handle to editPicture (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% global Project SYSINDEX Model_dir
% Project.System.TES(SYSINDEX).Picture=get(hObject,'string');
% axes(findobj(gcf,'tag','axesImage'));
% image(imread([Model_dir filesep 'graphics' filesep Project.System.TES(SYSINDEX).Picture]))
% set(gca,'tag','axesImage')
% axis off

% --- Executes during object creation, after setting all properties.
function editPicture_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editPicture (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbuttonPicture.
function pushbuttonPicture_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonPicture (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Project SYSINDEX 
[f,~]=uigetfile(Project.System.TES(SYSINDEX).Picture,'Select Picture File',which(Project.System.TES(SYSINDEX).Picture));
if f==0; return; end
Project.System.TES(SYSINDEX).Picture=f;
set(handles.editPicture,'string',Project.System.TES(SYSINDEX).Picture)
% editPicture_Callback(handles.editPicture,eventdata,handles)

% --- Executes on button press in pushbuttonSaveOnlyToProject.
function pushbuttonSaveOnlyToProject_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonSaveOnlyToProject (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clear global SYSINDEX TESoriginal
close(gcf)
uiresume

% --- Executes on button press in pushbuttonCancel.
function pushbuttonCancel_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonCancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Project TESoriginal
Project.System.TES = TESoriginal;
clear global SYSINDEX TESoriginal
close(gcf)

% --- Executes on button press in pushbuttonSaveAs.
function pushbuttonSaveAs_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonSaveAs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Project SYSINDEX Model_dir
savedir=fullfile(Model_dir,'System Library','ThermalStorage');
[f,p]=uiputfile(savedir,'Save As Thermal Storage Component',Project.System.TES(SYSINDEX).Name);
if f==0;return;end
component=Project.System.TES(SYSINDEX);
save([p f],'component')
clear global SYSINDEX TESoriginal
close(gcf)
uiresume

% --- Executes during object creation, after setting all properties.
function uipanelTEStype_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uipanelTEStype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



% --- Executes when selected object is changed in uipanelTEStype.
function uipanelTEStype_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanelTEStype 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
global Project SYSINDEX
switch get(eventdata.NewValue,'Tag')
    case 'radiobuttonColdTES'
        Project.System.TES(SYSINDEX).TEStype = 'cold';
    case 'radiobuttonHotTES'
        Project.System.TES(SYSINDEX).TEStype = 'hot';
end


function editSizeHours_Callback(hObject, eventdata, handles)
% hObject    handle to editSizeHours (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Project SYSINDEX
Project.System.TES(SYSINDEX).SizeHours = str2double(get(hObject,'string'));
OptimalTESsizing('editSizeHours');
COMPONENT = Project.System.TES(SYSINDEX);
set(handles.editSizeGal,'string',COMPONENT.SizeGal)
set(handles.editSize,'string',COMPONENT.SysSize)
set(handles.editOptSize,'string',COMPONENT.OptSize)
set(handles.editSizeHours,'string',COMPONENT.SizeHours)
fillRatesStatic(COMPONENT,handles)


% --- Executes during object creation, after setting all properties.
function editSizeHours_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editSizeHours (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editSizeGal_Callback(hObject, eventdata, handles)
% hObject    handle to editSizeGal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Project SYSINDEX
Project.System.TES(SYSINDEX).SizeGal = str2double(get(hObject,'string'));
Tcold = Project.System.TES(SYSINDEX).Tcold;
Thot = Project.System.TES(SYSINDEX).Thot;
kWh2GalDeg = 3600/4.186*.264; %3600 seconds / 4.186 kJ/kg*K * .264 gal/L
Project.System.TES(SYSINDEX).SysSize = Project.System.TES(SYSINDEX).SizeGal/kWh2GalDeg*(Thot-Tcold);
set(handles.editSize,'string',Project.System.TES(SYSINDEX).SysSize)
editSize_Callback(handles.editSize, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function editSizeGal_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editSizeGal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editTcold_Callback(hObject, eventdata, handles)
% hObject    handle to editTcold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Project SYSINDEX
Project.System.TES(SYSINDEX).Tcold = str2double(get(hObject,'string'));
editSize_Callback(handles.editSize, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function editTcold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editTcold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function editThot_Callback(hObject, eventdata, handles)
% hObject    handle to editThot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Project SYSINDEX
Project.System.TES(SYSINDEX).Thot = str2double(get(hObject,'string'));
editSize_Callback(handles.editSize, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function editThot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editThot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function editFillRate_Callback(hObject, eventdata, handles)
% hObject    handle to editFillRate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Project SYSINDEX
Project.System.TES(SYSINDEX).FillRate = str2double(get(hObject,'string'));
Project.System.TES(SYSINDEX).FillRatePerc = Project.System.TES(SYSINDEX).FillRate*60/Project.System.TES(SYSINDEX).SizeGal; % Fill rate in gal/min
set(handles.editFillRatePerc,'string',Project.System.TES(SYSINDEX).FillRatePerc)

% --- Executes during object creation, after setting all properties.
function editFillRate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editFillRate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function editFillRatePerc_Callback(hObject, eventdata, handles)
% hObject    handle to editFillRatePerc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Project SYSINDEX
Project.System.TES(SYSINDEX).FillRatePerc = str2double(get(hObject,'string'));
Project.System.TES(SYSINDEX).FillRate = Project.System.TES(SYSINDEX).SizeGal*Project.System.TES(SYSINDEX).FillRatePerc/60; % Fill rate in gal/min
set(handles.editFillRate,'string',Project.System.TES(SYSINDEX).FillRate)

% --- Executes during object creation, after setting all properties.
function editFillRatePerc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editFillRatePerc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editDischRate_Callback(hObject, eventdata, handles)
% hObject    handle to editDischRate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Project SYSINDEX
Project.System.TES(SYSINDEX).DischRate = str2double(get(hObject,'string'));
Project.System.TES(SYSINDEX).DischRatePerc = Project.System.TES(SYSINDEX).DischRate*60/Project.System.TES(SYSINDEX).SizeGal; % disch rate in gal/min
set(handles.editDischRatePerc,'string',Project.System.TES(SYSINDEX).DischRatePerc)

% --- Executes during object creation, after setting all properties.
function editDischRate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editDischRate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function editDischRatePerc_Callback(hObject, eventdata, handles)
% hObject    handle to editDischRatePerc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Project SYSINDEX
Project.System.TES(SYSINDEX).DischRatePerc = str2double(get(hObject,'string'));
Project.System.TES(SYSINDEX).DischRate = Project.System.TES(SYSINDEX).SizeGal*Project.System.TES(SYSINDEX).DischRatePerc/60; % disch rate in gal/min
set(handles.editDischRate,'string',Project.System.TES(SYSINDEX).DischRate)


% --- Executes during object creation, after setting all properties.
function editDischRatePerc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editDischRatePerc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function editLossDay_Callback(hObject, eventdata, handles)
% hObject    handle to editLossDay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Project SYSINDEX
Project.System.TES(SYSINDEX).LossDay = str2double(get(hObject,'string'));

% --- Executes during object creation, after setting all properties.
function editLossDay_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editLossDay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function fillRatesStatic(COMPONENT,handles)
global Project
demand = Project.Building;

% calculate the TES fill and discharge rates using demand & chiller capacity
FillRatePerc = 100/COMPONENT.SizeHours; % 100% divided by # of hours to fill
FillRate = COMPONENT.SizeGal*FillRatePerc/60; % fill rate in gal/min
if strcmp(COMPONENT.TEStype, 'cold')
    MaxDemandC = max(demand.DemandC);
    DischRatePerc = MaxDemandC/COMPONENT.SysSize*100;
    DischRate = COMPONENT.SizeGal*DischRatePerc/60; % disch rate in gal/min
elseif strcmp(COMPONENT.TEStype, 'hot')
    MaxDemandH = max(demand.DemandH); 
    DischRatePerc = MaxDemandH/COMPONENT.SysSize*100;
    DischRate = COMPONENT.SizeGal*DischRatePerc/60; % disch rate in gal/min
end
set(handles.textFillRate,'string',FillRate)
set(handles.textFillRatePerc,'string',FillRatePerc)
set(handles.textDischRate,'string',DischRate)
set(handles.textDischRatePerc,'string',DischRatePerc)
