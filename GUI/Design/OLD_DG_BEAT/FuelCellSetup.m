function varargout = FuelCellSetup(varargin)
% FUELCELLSETUP MATLAB code for FuelCellSetup.fig
%      FUELCELLSETUP, by itself, creates a new FUELCELLSETUP or raises the existing
%      singleton*.
%
%      H = FUELCELLSETUP returns the handle to a new FUELCELLSETUP or the handle to
%      the existing singleton*.
%
%      FUELCELLSETUP('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FUELCELLSETUP.M with the given input arguments.
%
%      FUELCELLSETUP('Property','Value',...) creates a new FUELCELLSETUP or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before FuelCellSetup_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to FuelCellSetup_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help FuelCellSetup

% Last Modified by GUIDE v2.5 26-Mar-2013 12:05:32

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @FuelCellSetup_OpeningFcn, ...
                   'gui_OutputFcn',  @FuelCellSetup_OutputFcn, ...
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


% --- Executes just before FuelCellSetup is made visible.
function FuelCellSetup_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to FuelCellSetup (see VARARGIN)

% Add a wait dialogue box
figure_FC = get(0,'CurrentFigure');
MSG_FC=msgbox('Loading');
figure(figure_FC)

% Choose default command line output for FuelCellSetup
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes FuelCellSetup wait for user response (see UIRESUME)
% uiwait(handles.figure1);
global Project SYSINDEX StackYrs CHPoriginal StackLifekWh
StackLifekWh =Project.Economic.LifekWh;
CHPoriginal = Project.System.CHP(SYSINDEX);
OptimalFCsizing(SYSINDEX,'editMaxPower');
COMPONENT=Project.System.CHP(SYSINDEX);
StackYrs = Project.Economic.LifeYrs;
COMPONENT.TurnDown = COMPONENT.SysSize(1,1)/COMPONENT.SysSize(2,1);

set(handles.editName,'string',COMPONENT.Name)
set(handles.textType,'string',COMPONENT.Type)
set(handles.editVersion,'string',COMPONENT.Version)
set(handles.editDescription,'string',COMPONENT.Description)
% set(handles.editPicture,'string',COMPONENT.Picture)
% editPicture_Callback(handles.editPicture,eventdata,handles)

set(handles.editMinPower,'string',round(COMPONENT.SysSize(2,1)))
set(handles.editMaxPower,'string',round(COMPONENT.SysSize(1,1)))
set(handles.editResponseTime,'string',COMPONENT.SysRamp)
set(handles.editTurnDown,'string',COMPONENT.TurnDown)
set(handles.OptSize,'string',round(COMPONENT.OptSize))
set(handles.uitableCHP,'Data',COMPONENT.SysCHP)
set(handles.uitableEfficiency,'Data',COMPONENT.SysEff)
updatePlot()

close(MSG_FC)

% --- Outputs from this function are returned to the command line.
function varargout = FuelCellSetup_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes when entered data in editable cell(s) in uitableEfficiency.
function uitableEfficiency_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to uitableEfficiency (see GCBO)
% eventdata  structure with the following fields (see UITABLEEFFICIENCY)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
global Project SYSINDEX
Project.System.CHP(SYSINDEX).SysEff = get(handles.uitableEfficiency,'Data');
updatePlot()

% --- Executes when selected cell(s) is changed in uitableEfficiency.
function uitableEfficiency_CellSelectionCallback(hObject, eventdata, handles)
% hObject    handle to uitableEfficiency (see GCBO)
% eventdata  structure with the following fields (see UITABLEEFFICIENCY)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)

function editMinPower_Callback(hObject, eventdata, handles)
% hObject    handle to editMinPower (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Project SYSINDEX
Project.System.CHP(SYSINDEX).SysSize(2,1)=str2double(get(hObject,'string'));
Project.System.CHP(SYSINDEX).TurnDown = Project.System.CHP(SYSINDEX).SysSize(1,1)/Project.System.CHP(SYSINDEX).SysSize(2,1);
set(handles.editTurnDown,'string',Project.System.CHP(SYSINDEX).TurnDown)


% --- Executes during object creation, after setting all properties.
function editMinPower_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editMinPower (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editResponseTime_Callback(hObject, eventdata, handles)
% hObject    handle to editResponseTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Project SYSINDEX
Project.System.CHP(SYSINDEX).SysRamp=str2double(get(hObject,'string'));


% --- Executes during object creation, after setting all properties.
function editResponseTime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editResponseTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pushbuttonSaveOnlyToProject.
function pushbuttonSaveOnlyToProject_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonSaveOnlyToProject (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Project StackLifekWh
Project.Economic.LifekWh  = StackLifekWh;
clear global SYSINDEX CHPoriginal StackLifekWh StackYrs
close(gcf)
uiresume

% --- Executes on button press in pushbuttonCancel.
function pushbuttonCancel_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonCancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Project SYSINDEX CHPoriginal
Project.System.CHP(SYSINDEX) = CHPoriginal;
clear global SYSINDEX CHPoriginal StackLifekWh StackYrs
close(gcf)


function editName_Callback(hObject, eventdata, handles)
% hObject    handle to editName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Project SYSINDEX
Project.System.CHP(SYSINDEX).Name=get(hObject,'string');

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
Project.System.CHP(SYSINDEX).Description=get(hObject,'string');

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
Project.System.CHP(SYSINDEX).Version=get(hObject,'string');

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
% Project.System.CHP(SYSINDEX).Picture=get(hObject,'string');
% axes(findobj(gcf,'tag','axesImage'));
% image(imread([Model_dir filesep 'graphics' filesep Project.System.CHP(SYSINDEX).Picture]))
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
f=uigetfile(Project.System.CHP(SYSINDEX).Picture,'Select Picture File',which(Project.System.CHP(SYSINDEX).Picture));
if f==0; return; end
Project.System.CHP(SYSINDEX).Picture=f;
set(handles.editPicture,'string',Project.System.CHP(SYSINDEX).Picture)
% editPicture_Callback(handles.editPicture,eventdata,handles)

% --- Executes on button press in pushbuttonSaveAs.
function pushbuttonSaveAs_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonSaveAs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Project SYSINDEX  StackLifekWh Model_dir
[f,p]=uiputfile('*.mat','Save As CHP Component',fullfile(Model_dir,'System Library','CHP'));
if f==0;return;end
component=Project.System.CHP(SYSINDEX);
save([p f],'component')
Project.Economic.LifekWh  = StackLifekWh;
clear global SYSINDEX CHPoriginal StackLifekWh StackYrs
close(gcf)
uiresume


function editMaxPower_Callback(hObject, eventdata, handles)
% hObject    handle to editMaxPower (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Project SYSINDEX StackYrs StackLifekWh
Project.System.CHP(SYSINDEX).SysSize(1,1)=str2double(get(hObject,'string'));
OptimalFCsizing(SYSINDEX,'editMaxPower');
set(handles.editMaxPower,'string',round(Project.System.CHP(SYSINDEX).SysSize(1,1)))
set(handles.OptSize,'string',round(Project.System.CHP(SYSINDEX).OptSize))
StackLifekWh = StackYrs*Project.System.CHP(SYSINDEX).SysSize(1,1)*8760;

% --- Executes during object creation, after setting all properties.
function editMaxPower_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editMaxPower (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when entered data in editable cell(s) in uitableCHP.
function uitableCHP_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to uitableCHP (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
global Project SYSINDEX 
Project.System.CHP(SYSINDEX).SysCHP=get(handles.uitableCHP,'Data');
axes(findobj(gcf,'tag','axes2'));
plot(Project.System.CHP(SYSINDEX).SysCHP(:,1),Project.System.CHP(SYSINDEX).SysCHP(:,2),'b-o');
set(gca,'tag','axes2')
ylim([0 1])
xlabel('Temp [deg C]')
ylabel('Derate Factor')
title('System CHP')


function OptSize_Callback(hObject, eventdata, handles)
% hObject    handle to OptSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Project SYSINDEX StackYrs StackLifekWh
Project.System.CHP(SYSINDEX).OptSize=str2double(get(hObject,'string'));
OptimalFCsizing(SYSINDEX,'OptSize');
set(handles.editMaxPower,'string',round(Project.System.CHP(SYSINDEX).SysSize(1,1)))
set(handles.editMinPower,'string',round(Project.System.CHP(SYSINDEX).SysSize(2,1)))
set(handles.OptSize,'string',round(Project.System.CHP(SYSINDEX).OptSize))
StackLifekWh = StackYrs*Project.System.CHP(SYSINDEX).SysSize(1,1)*8760;

% --- Executes during object creation, after setting all properties.
function OptSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to OptSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function editTurnDown_Callback(hObject, eventdata, handles)
% hObject    handle to editTurnDown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Project SYSINDEX 
Project.System.CHP(SYSINDEX).TurnDown=str2double(get(hObject,'string'));
Project.System.CHP(SYSINDEX).SysSize(2,1) = Project.System.CHP(SYSINDEX).SysSize(1,1)/Project.System.CHP(SYSINDEX).TurnDown;
set(handles.editMinPower,'string',round(Project.System.CHP(SYSINDEX).SysSize(2,1)))

% --- Executes during object creation, after setting all properties.
function editTurnDown_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editTurnDown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function updatePlot()
global Project SYSINDEX
X = Project.System.CHP(SYSINDEX).SysEff(:,1);
axes(findobj(gcf,'tag','axes1'))
[Ax,H1,H2] = plotyy(X,Project.System.CHP(SYSINDEX).SysEff(:,2),X,Project.System.CHP(SYSINDEX).SysEff(:,4));
set(H1, 'LineWidth',3,'Color','blue','LineStyle','-')
set(H2, 'LineWidth',3,'Color','green','LineStyle','--')
set(Ax(1),'XColor','k','YColor','k','Fontsize',12)
set(Ax(2),'XColor','k','YColor','k','Fontsize',12)
hold(Ax(1));
hold(Ax(2));
plot(Ax(1), X,Project.System.CHP(SYSINDEX).SysEff(:,3),'r-o','LineWidth',2)
plot(Ax(2), X,Project.System.CHP(SYSINDEX).SysEff(:,5),'c-.','LineWidth',2)
plot(Ax(2), X,Project.System.CHP(SYSINDEX).SysEff(:,6),'m--','LineWidth',2)

set(gca,'tag','axes1')
axis(Ax(1),[0 1 0 1])
axis(Ax(2),[0 1 0 2000])
set(Ax(1),'YTick',0:.1:1)
set(Ax(2),'YTick',0:200:2000)
ylabel(Ax(1),'Efficiency [%]','Fontsize',16,'Color',[0 0 0])
ylabel(Ax(2),'Emissions (lb/MWh)','Fontsize',16,'Color',[0 0 0])
legend(Ax(1),'Electric Efficiency (% LHV)', 'Heat Recovery (% of fuel energy)','CO2 emissions','NOx emissions','SO2 emissions','Location','Best')
