function varargout = ChillerSetup(varargin)
% CHILLERSETUP M-file for ChillerSetup.fig
%      CHILLERSETUP, by itself, creates a new CHILLERSETUP or raises the existing
%      singleton*.
%
%      H = CHILLERSETUP returns the handle to a new CHILLERSETUP or the handle to
%      the existing singleton*.
%
%      CHILLERSETUP('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CHILLERSETUP.M with the given input arguments.
%
%      CHILLERSETUP('Property','Value',...) creates a new CHILLERSETUP or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ChillerSetup_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ChillerSetup_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ChillerSetup

% Last Modified by GUIDE v2.5 17-Jan-2013 13:36:12

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ChillerSetup_OpeningFcn, ...
                   'gui_OutputFcn',  @ChillerSetup_OutputFcn, ...
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


% --- Executes just before ChillerSetup is made visible.
function ChillerSetup_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ChillerSetup (see VARARGIN)

% Add a wait dialogue box
figure_Chill = get(0,'CurrentFigure');
MSG_Chill=msgbox('Loading');
figure(figure_Chill)

% Choose default command line output for ChillerSetup
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

global Project SYSINDEX ChillOriginal
ChillOriginal = Project.System.Chiller;
OptimalChillerSizing('editSize');
kWh2ton = 1/3.51685;
COMPONENT=Project.System.Chiller(SYSINDEX);
COMPONENT.SizeTons = COMPONENT.SysSize*kWh2ton;

set(handles.editName,'string',COMPONENT.Name)
set(handles.textType,'string',COMPONENT.Type)
set(handles.editVersion,'string',COMPONENT.Version)
set(handles.editDescription,'string',COMPONENT.Description)
% set(handles.editPicture,'string',COMPONENT.Picture)
% editPicture_Callback(handles.editPicture,eventdata,handles)
set(handles.editSize,'string',COMPONENT.SysSize)
set(handles.editSizeTons,'string',COMPONENT.SizeTons)
set(handles.editCOP,'string',COMPONENT.COP)
set(handles.editOptSize,'string',COMPONENT.OptSize)
set(handles.uitableCOPcurve,'Data',COMPONENT.COPcurve);
ButtonMenu = handles.uipanelChillType;
radio1 = handles.radiobuttonElecChill;
radio2 = handles.radiobuttonAbsorpChill;
if strcmp(COMPONENT.ChillType, 'electric')
    set(ButtonMenu,'SelectedObject', radio1)
    set(handles.textOptSize,'string','% Max Cool Demand')
else
    set(ButtonMenu,'SelectedObject', radio2)
    set(handles.textOptSize,'string','% Max Heat Available')
end

axes(findobj(gcf,'tag','axesCOPplot'));
h(1)=plot(COMPONENT.COPcurve(:,1),(COMPONENT.COPcurve(:,2))*COMPONENT.COP,'b-o');
set(gca,'tag','axesCOPplot')
ylim([0 COMPONENT.COP])
xlabel('Power [%]')
ylabel('Chiller Coefficient of Performance')
title('Chiller Efficiency')

close(MSG_Chill)


% --- Outputs from this function are returned to the command line.
function varargout = ChillerSetup_OutputFcn(hObject, eventdata, handles) 
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
Project.System.Chiller(SYSINDEX).SysSize=str2double(get(hObject,'string'));
OptimalChillerSizing('editSize');
kWh2ton = 1/3.51685;
set(handles.editOptSize,'string',Project.System.Chiller(SYSINDEX).OptSize)
Project.System.Chiller(SYSINDEX).SizeTons = Project.System.Chiller(SYSINDEX).SysSize*kWh2ton;
set(handles.editSizeTons,'string',Project.System.Chiller(SYSINDEX).SizeTons)


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



function editCOP_Callback(hObject, eventdata, handles)
% hObject    handle to editCOP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editCOP as text
%        str2double(get(hObject,'String')) returns contents of editCOP as a double
global Project SYSINDEX
Project.System.Chiller(SYSINDEX).COP=str2double(get(hObject,'string'));
COMPONENT=Project.System.Chiller(SYSINDEX);
axes(findobj(gcf,'tag','axesCOPplot'));
h(1)=plot(COMPONENT.COPcurve(:,1),(COMPONENT.COPcurve(:,2))*COMPONENT.COP,'b-o');
set(gca,'tag','axesCOPplot')
ylim([0 COMPONENT.COP])
xlabel('Power [%]')
ylabel('Chiller Coefficient of Performance')
title('Chiller Efficiency')

% --- Executes during object creation, after setting all properties.
function editCOP_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editCOP (see GCBO)
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
Project.System.Chiller(SYSINDEX).OptSize=str2double(get(hObject,'string'));
OptimalChillerSizing('OptSize');
kWh2ton = 1/3.51685;
set(handles.editSize,'string',Project.System.Chiller(SYSINDEX).SysSize)
Project.System.Chiller(SYSINDEX).SizeTons = Project.System.Chiller(SYSINDEX).SysSize*kWh2ton;
set(handles.editSizeTons,'string',Project.System.Chiller(SYSINDEX).SizeTons)

% --- Executes during object creation, after setting all properties.
function editOptSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editOptSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function editName_Callback(hObject, eventdata, handles)
% hObject    handle to editName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editName as text
%        str2double(get(hObject,'String')) returns contents of editName as a double
global Project SYSINDEX
Project.System.Chiller(SYSINDEX).Name=get(hObject,'string');

% --- Executes during object creation, after setting all properties.
function editName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editDescription_Callback(hObject, eventdata, handles)
% hObject    handle to editDescription (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editDescription as text
%        str2double(get(hObject,'String')) returns contents of editDescription as a double
global Project SYSINDEX
Project.System.Chiller(SYSINDEX).Description=get(hObject,'string');

% --- Executes during object creation, after setting all properties.
function editDescription_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editDescription (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editVersion_Callback(hObject, eventdata, handles)
% hObject    handle to editVersion (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editVersion as text
%        str2double(get(hObject,'String')) returns contents of editVersion as a double
global Project SYSINDEX
Project.System.Chiller(SYSINDEX).Version=get(hObject,'string');

% --- Executes during object creation, after setting all properties.
function editVersion_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editVersion (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editPicture_Callback(hObject, eventdata, handles)
% hObject    handle to editPicture (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editPicture as text
%        str2double(get(hObject,'String')) returns contents of editPicture as a double
% global Project SYSINDEX Model_dir
% Project.System.Chiller(SYSINDEX).Picture=get(hObject,'string');
% axes(findobj(gcf,'tag','axesImage'));
% image(imread([Model_dir filesep 'graphics' filesep Project.System.Chiller(SYSINDEX).Picture]))
% set(gca,'tag','axesImage')
% axis off

% --- Executes during object creation, after setting all properties.
function editPicture_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editPicture (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbuttonPicture.
function pushbuttonPicture_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonPicture (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Project SYSINDEX
COMPONENT = Project.System.Chiller(SYSINDEX);
[f,~]=uigetfile(COMPONENT.Picture,'Select Picture File',which(COMPONENT.Picture));
if f==0; return; end
COMPONENT.Picture=f;
set(handles.editPicture,'string',COMPONENT.Picture)
% editPicture_Callback(handles.editPicture,eventdata,handles)

% --- Executes on button press in pushbuttonSaveOnlyToProject.
function pushbuttonSaveOnlyToProject_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonSaveOnlyToProject (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clear global SYSINDEX ChillOriginal
close(gcf)
uiresume

% --- Executes on button press in pushbuttonCancel.
function pushbuttonCancel_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonCancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Project ChillOriginal
Project.System.Chiller = ChillOriginal;
clear global SYSINDEX ChillOriginal
close(gcf)

% --- Executes on button press in pushbuttonSaveAs.
function pushbuttonSaveAs_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonSaveAs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Project SYSINDEX Model_dir
[f,p]=uiputfile('*.mat','Save As Chiller Component',fullfile(Model_dir,'System Library','Chiller'));
if f==0;return;end
component=Project.System.Chiller(SYSINDEX);
save([p f],'component')
clear global SYSINDEX ChillOriginal
close(gcf)
uiresume

% --- Executes during object creation, after setting all properties.
function uipanelChillType_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uipanelChillType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes when entered data in editable cell(s) in uitableCOPcurve.
function uitableCOPcurve_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to uitableCOPcurve (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
global Project SYSINDEX
Data=get(handles.uitableCOPcurve,'Data');
Project.System.Chiller(SYSINDEX).COPcurve=Data;
COMPONENT = Project.System.Chiller(SYSINDEX);
set(handles.uitableCOPcurve,'Data',COMPONENT.COPcurve);
axes(findobj(gcf,'tag','axesCOPplot'))
h(1)=plot(COMPONENT.COPcurve(:,1),COMPONENT.COPcurve(:,2)*COMPONENT.COP,'b-o');
set(gca,'tag','axesCOPplot')
ylim([0 COMPONENT.COP])
xlabel('Power [%]')
ylabel('Chiller Coefficient of Performance')
title('Chiller Efficiency')


% --- Executes when selected object is changed in uipanelChillType.
function uipanelChillType_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanelChillType 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
global Project SYSINDEX
switch get(eventdata.NewValue,'Tag')
    case 'radiobuttonElecChill'
        Project.System.Chiller(SYSINDEX).ChillType = 'electric';
        set(handles.textOptSize,'string','% Max Cool Demand')
    case 'radiobuttonAbsorpChill'
        Project.System.Chiller(SYSINDEX).ChillType = 'absorption';
        set(handles.textOptSize,'string','% Max Heat Available')
    otherwise
        Project.System.Chiller(SYSINDEX).ChillType = 'electric';
end


% --- Executes when selected cell(s) is changed in uitableCOPcurve.
function uitableCOPcurve_CellSelectionCallback(hObject, eventdata, handles)
% hObject    handle to uitableCOPcurve (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function uitableCOPcurve_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uitableCOPcurve (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


function editSizeTons_Callback(hObject, eventdata, handles)
% hObject    handle to editSizeTons (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editSizeTons as text
%        str2double(get(hObject,'String')) returns contents of editSizeTons as a double
global Project SYSINDEX
Project.System.Chiller(SYSINDEX).SizeTons = str2double(get(hObject,'string'));
kWh2ton = 1/3.51685;
Project.System.Chiller(SYSINDEX).SysSize = Project.System.Chiller(SYSINDEX).SizeTons/kWh2ton;
set(handles.editSize,'string',Project.System.Chiller(SYSINDEX).SysSize)
editSize_Callback(handles.editSize, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function editSizeTons_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editSizeTons (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
