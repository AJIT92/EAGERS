function varargout = SetupGenerator(varargin)
% SETUPGENERATOR MATLAB code for SetupGenerator.fig
%      SETUPGENERATOR, by itself, creates a new SETUPGENERATOR or raises the existing
%      singleton*.
%
%      H = SETUPGENERATOR returns the handle to a new SETUPGENERATOR or the handle to
%      the existing singleton*.
%
%      SETUPGENERATOR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SETUPGENERATOR.M with the given input arguments.
%
%      SETUPGENERATOR('Property','Value',...) creates a new SETUPGENERATOR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SetupGenerator_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SetupGenerator_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SetupGenerator

% Last Modified by GUIDE v2.5 23-May-2015 16:25:31

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SetupGenerator_OpeningFcn, ...
                   'gui_OutputFcn',  @SetupGenerator_OutputFcn, ...
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

% --- Executes just before SetupGenerator is made visible.
function SetupGenerator_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SetupGenerator (see VARARGIN)

% Choose default command line output for SetupGenerator
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

global Plant SYSINDEX Generator figHandle
figHandle = handles.figure1;

Generator = Plant.Generator(SYSINDEX);

set(handles.editName,'string',Generator.Name)
if max(Generator.Output.Cooling)>0
    set(handles.editSize,'string',Generator.Size./3.517)
    set(handles.editSize,'Value',Generator.Size./3.517)
else
    set(handles.editSize,'string',Generator.Size)
end

if strcmp(Generator.Source, 'NG')
    set(handles.uipanelGenType,'SelectedObject', handles.radiobuttonNatGas)
elseif strcmp(Generator.Source, 'Oil')
    set(handles.uipanelGenType,'SelectedObject', handles.radiobuttonOil)
elseif strcmp(Generator.Source, 'Electricity')
    set(handles.uipanelGenType,'SelectedObject', handles.radiobuttonElectricity)
elseif strcmp(Generator.Source, 'Heat')
    set(handles.uipanelGenType,'SelectedObject', handles.radiobuttonHeat)
elseif strcmp(Generator.Source, 'Steam')
    set(handles.uipanelGenType,'SelectedObject', handles.radiobuttonSteam)
elseif strcmp(Generator.Source, 'Biofuel')
    set(handles.uipanelGenType,'SelectedObject', handles.radiobuttonBiofuel)
end
if max(Generator.Output.Electricity)>0
    set(handles.checkboxElec,'Value',1)
end
if max(Generator.Output.Steam)>0
    set(handles.checkboxSteam,'Value',1)
end
if max(Generator.Output.Heat)>0
    set(handles.checkboxHeat,'Value',1)
end
if max(Generator.Output.Cooling)>0
    set(handles.checkboxCooling,'Value',1)
end
setSizeUnits(handles)
plotGenEfficiency(Generator,handles)
plotStartup(handles)
A = Generator.VariableStruct.StateSpace.A;
B = Generator.VariableStruct.StateSpace.B;
C = Generator.VariableStruct.StateSpace.C;
D = Generator.VariableStruct.StateSpace.D;
if isfield(Generator.VariableStruct.StateSpace,'Dt')
    Dt = Generator.VariableStruct.StateSpace.Dt;
else Dt =1;
    Generator.VariableStruct.StateSpace.Dt = 1;
end
SS = ss(A,B,C,D,Dt);
[Wn,zeta] = damp(SS);
set(handles.editFreq,'string',num2str(Wn(1)/(2*pi()/1000)))
set(handles.editDamp,'string',num2str(zeta(1)))

if isfield(Generator.VariableStruct,'RestartTime')
    RestartTime = Generator.VariableStruct.RestartTime;
else RestartTime =0;
end
set(handles.editRestartTime,'string',num2str(RestartTime));
plotResponse(handles,A,B,C,D,Dt,'Discrete')

% --- Outputs from this function are returned to the command line.
function varargout = SetupGenerator_OutputFcn(hObject, eventdata, handles) 
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
global Generator
Generator.Name=get(hObject,'string');

function editName_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function editSize_Callback(hObject, eventdata, handles)
% hObject    handle to editSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Generator
Generator.Size=str2double(get(hObject,'string'));
set(handles.editSize,'Value',str2double(get(hObject,'string')));


function editSize_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function setSizeUnits(handles)
if get(handles.checkboxElec,'Value')==1
    set(handles.textUnits,'string','kW')
elseif get(handles.checkboxCooling,'Value')==1
    set(handles.textUnits,'string','tons')
elseif get(handles.checkboxSteam,'Value')==1
    set(handles.textUnits,'string','pph')
elseif get(handles.checkboxHeat,'Value')==1
    set(handles.textUnits,'string','kW')
end
%count the types of outputs
Products = get(handles.checkboxElec,'Value')+get(handles.checkboxCooling,'Value')+get(handles.checkboxSteam,'Value')+get(handles.checkboxHeat,'Value');
if Products >1
    set(handles.editFreq,'Enable','off')
    set(handles.editDamp,'Enable','off')
else 
    set(handles.editFreq,'Enable','on')
    set(handles.editDamp,'Enable','on')
end


% --- Executes on button press in checkboxElec.
function checkboxElec_Callback(hObject, eventdata, handles)
setSizeUnits(handles)
pushbuttonEditTable_Callback(hObject, eventdata, handles)

% --- Executes on button press in checkboxSteam.
function checkboxSteam_Callback(hObject, eventdata, handles)
setSizeUnits(handles)
pushbuttonEditTable_Callback(hObject, eventdata, handles)
% --- Executes on button press in checkboxHeat.
function checkboxHeat_Callback(hObject, eventdata, handles)
setSizeUnits(handles)
pushbuttonEditTable_Callback(hObject, eventdata, handles)
% --- Executes on button press in checkboxCooling.
function checkboxCooling_Callback(hObject, eventdata, handles)
setSizeUnits(handles)
pushbuttonEditTable_Callback(hObject, eventdata, handles)
% --- Executes during object creation, after setting all properties.
function uipanelGenType_CreateFcn(hObject, eventdata, handles)

% --- Executes when selected object is changed in uipanelGenType.
function uipanelGenType_SelectionChangeFcn(hObject, eventdata, handles)
global Generator
switch get(eventdata.NewValue,'Tag')
    case 'radiobuttonNatGas'
        Generator.Source = 'NG';
    case 'radiobuttonElectricity'
        Generator.Source = 'Electricity';
    case 'radiobuttonHeat'
        Generator.Source = 'Heat';
    case 'radiobuttonSteam'
        Generator.Source = 'Steam';
    case 'radiobuttonOil'
        Generator.Source = 'Oil';
    case 'radiobuttonBiofuel'
        Generator.Source = 'Biofuel';
end


% --- Executes on button press in pushbuttonEditTable.
function pushbuttonEditTable_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonEditTable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Generator EditTableHandle TableData
ColNames = {'Capacity'};
TableData =Generator.Output.Capacity;
if get(handles.checkboxElec,'Value')==1
    ColNames(end+1) = {'Electricity'};
    TableData(:,end+1) = Generator.Output.Electricity;
else Generator.Output.Electricity = zeros(length(Generator.Output.Capacity),1);
end
if get(handles.checkboxSteam,'Value')==1
    ColNames(end+1) = {'Steam'};
    TableData(:,end+1) = Generator.Output.Steam;
else Generator.Output.Steam = zeros(length(Generator.Output.Capacity),1);
end
if get(handles.checkboxHeat,'Value')==1
    ColNames(end+1) = {'Heat'};
    TableData(:,end+1) = Generator.Output.Heat;
else Generator.Output.Heat = zeros(length(Generator.Output.Capacity),1);
end
if get(handles.checkboxCooling,'Value')==1
    ColNames(end+1) = {'Cooling'};
    TableData(:,end+1) = Generator.Output.Cooling;
else Generator.Output.Cooling = zeros(length(Generator.Output.Capacity),1);
end
Text1 = 'Values represent output energy product per unit of energy source input';
EditTable(ColNames,Text1);
waitfor(EditTableHandle)
Data = TableData;
Generator.Output.Capacity = Data(:,1);
if get(handles.checkboxElec,'Value')==1
    Generator.Output.Electricity = Data(:,2);
else Generator.Output.Electricity = 0*Data(:,1);
end
if get(handles.checkboxSteam,'Value')==1
    c = find(strcmp(ColNames,'Steam'));
    Generator.Output.Steam = Data(:,c);
else Generator.Output.Steam = 0*Data(:,1);
end
if get(handles.checkboxHeat,'Value')==1
    c = find(strcmp(ColNames,'Heat'));
    Generator.Output.Heat = Data(:,c);
else Generator.Output.Heat = 0*Data(:,1);
end
if get(handles.checkboxCooling,'Value')==1
    c = find(strcmp(ColNames,'Cooling'));
    Generator.Output.Cooling = Data(:,c);
else Generator.Output.Cooling = 0*Data(:,1);
end
plotGenEfficiency(Generator,handles)

function plotGenEfficiency(Gen,handles)
axes(handles.axes1)
hold off
str = {};
if get(handles.checkboxElec,'Value')==1
    c = Gen.Output.Capacity./Gen.Output.Electricity;
    [AX, H1, H2] = plotyy(Gen.Output.Capacity,Gen.Output.Electricity,Gen.Output.Capacity(2:end),c(2:end));
    hold on
    set(H1,'Color','k','LineStyle','-','LineWidth',2,'Marker','o')
    str(end+1) = {'Electric'};
end
if get(handles.checkboxHeat,'Value')==1
    if get(handles.checkboxElec,'Value')==1
       plot(AX(1),Gen.Output.Capacity,Gen.Output.Heat,'r-o');
    else c = Gen.Output.Capacity./Gen.Output.Heat; 
        [AX, H1, H2] = plotyy(Gen.Output.Capacity,Gen.Output.Heat,Gen.Output.Capacity(2:end),c(2:end));
        set(H1,'Color','r','LineStyle','-','LineWidth',2,'Marker','o')
        hold on
    end
    str(end+1) = {'Heat'};
end
if get(handles.checkboxSteam,'Value')==1
    c = Gen.Output.Capacity./Gen.Output.Steam; 
    [AX, H1, H2] = plotyy(Gen.Output.Capacity,Gen.Output.Steam,Gen.Output.Capacity(2:end),c(2:end));
    set(H1,'Color','m','LineStyle','-','LineWidth',2,'Marker','o')
    hold on
    str(end+1) = {'Steam'};
end
if get(handles.checkboxCooling,'Value')==1
    c = Gen.Output.Capacity./Gen.Output.Cooling; 
    [AX, H1, H2] = plotyy(Gen.Output.Capacity,Gen.Output.Cooling,Gen.Output.Capacity(2:end),c(2:end));
    set(H1,'Color','b','LineStyle','-','LineWidth',2,'Marker','o')
    str(end+1) = {'Cooling'};
end
set(AX,{'ycolor'},{'k';'k'})
set(H2,'Color','g','LineStyle',':','LineWidth',3)
str(end+1) = {'Cost Curve'};

set(gca,'tag','axes1')
if max(Gen.Output.Cooling)>1
    ylim([0 max(Gen.Output.Cooling)])
    a = round(max(Gen.Output.Cooling))+1;
    ylim(AX(1),[0,a])
    ylim(AX(2),[0,1])
else
    a = round(max(c))+1;
    ylim(AX(1),[0,1])
    ylim(AX(2),[0,a])
    set(AX(1),'YTick',0:.1:1)
    set(AX(2),'YTick',0:a/10:a)
end
set(get(AX(1),'Ylabel'),'String','Efficiency')
set(get(AX(2),'Ylabel'),'String','Cost Curve Shape')
xlabel('% of Capacity')
legend(str)
title('Efficiency / Cost')

% --- Executes on button press in pushbuttonSIMO.
%add back in when you can do single input multiple output response
function pushbuttonSIMO_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonSIMO (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global SIMOhandle Generator
SIMOresponse;
waitfor(SIMOhandle)
A = Generator.VariableStruct.StateSpace.A;
B = Generator.VariableStruct.StateSpace.B;
C = Generator.VariableStruct.StateSpace.C;
D = Generator.VariableStruct.StateSpace.D;
Dt = Generator.VariableStruct.StateSpace.Dt;%sampling time
SS = ss(A,B,C,D,Dt);
[Wn,zeta] = damp(SS);
set(handles.editFreq,'string',num2str(Wn(1)/(2*pi()/1000)))
set(handles.editDamp,'string',num2str(zeta(1)))
plotResponse(handles,A,B,C,D,Dt,'Discrete')

function editFreq_Callback(hObject, eventdata, handles)
w0 = str2double(get(handles.editFreq,'String'))*2*pi()/1000;
zeta = str2double(get(handles.editDamp,'String'));
A = [0 1; -(w0^2) -2*zeta*w0;];
B = [0; (w0^2);];
C = [1 0];
D = 0;
plotResponse(handles,A,B,C,D,[],'Continuous')

function editDamp_Callback(hObject, eventdata, handles)
w0 = str2double(get(handles.editFreq,'String'))*2*pi()/1000;
zeta = str2double(get(handles.editDamp,'String'));
A = [0 1; -(w0^2) -2*zeta*w0;];
B = [0; (w0^2);];
C = [1 0];
D = 0;
plotResponse(handles,A,B,C,D,[],'Continuous')

function plotResponse(handles, A, B, C, D,Dt,CorD)
global Generator
if strcmp(CorD,'Continuous')
    SS = ss(A,B,C,D);
    SS = c2d(SS,1);
    r = length(SS);
    A = SS(1).A;
    B = SS(1).B;
    C = SS(1).C;
    D = SS(1).D;
    for i = 2:1:r
        C(end+1,:) = SS(i).C;
        D(end+1,:) = SS(i).D;
    end
    Dt=1;
else SS = ss(A,B,C,D,Dt);
end
Generator.VariableStruct.StateSpace.A = A;
Generator.VariableStruct.StateSpace.B = B;
Generator.VariableStruct.StateSpace.C = C;
Generator.VariableStruct.StateSpace.D = D;
x0 = [];
Names = {};
if get(handles.checkboxElec,'Value')==1
    x0(end+1) = Generator.VariableStruct.Startup.Electricity(end);
    Names(end+1) = {'Electricity'};
end
if get(handles.checkboxCooling,'Value')==1
    x0(end+1) = Generator.VariableStruct.Startup.Cooling(end);
    Names(end+1) = {'Cooling'};
end
if get(handles.checkboxSteam,'Value')==1
    x0(end+1) = Generator.VariableStruct.Startup.Steam(end);
    Names(end+1) = {'Steam'};
end
if get(handles.checkboxHeat,'Value')==1
    x0(end+1) = Generator.VariableStruct.Startup.Heat(end);
    Names(end+1) = {'Heat'};
end
nS = round(3600/Dt)+1;
t = linspace(0, Dt*(nS-1),nS);
u = Generator.Size*linspace(1,1,nS);
[n,n2] = size(C);
X0 = zeros(n2,1);
for i = 1:1:n
    X0(find(C(i,:),1,'first'))=x0(i);
end
% r = size(SS(1,1).A);
% s = length(x0);
% x0 = [x0;zeros(r(1)-s,1);];
[y,t] = lsim(SS,u,t,X0);
plot(handles.axes2,t/60,y);
xlabel('Time (min)')
legend(Names)


function pushbuttonStartUp_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonStartUp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Generator EditTableHandle TableData
Start = Generator.VariableStruct.Startup;
TableData =[];
ColNames = {'Time'};
TableData(:,1) = Start.Time';
if get(handles.checkboxElec,'Value')==1
    ColNames(end+1) = {'Electricity'};
    TableData(:,end+1) = Start.Electricity';
else Start.Electricity = zeros(1,length(Start.Time));
end
if get(handles.checkboxSteam,'Value')==1
    ColNames(end+1) = {'Steam'};
    TableData(:,end+1) = Start.Steam';
else Start.Steam = zeros(1,length(Start.Time));
end
if get(handles.checkboxHeat,'Value')==1
    ColNames(end+1) = {'Heat'};
    TableData(:,end+1) = Start.Heat';
else Start.Heat = zeros(1,length(Start.Time));
end
if get(handles.checkboxCooling,'Value')==1
    ColNames(end+1) = {'Cooling'};
    TableData(:,end+1) = Start.Cooling';
else Start.Cooling = zeros(1,length(Start.Time));
end
ColNames(end+1) = {'Input'};
TableData(:,end+1) = Start.Input';
Text1 = 'Values represent normalized output 1 = 100%'; 
EditTable(ColNames,Text1);
waitfor(EditTableHandle)
Start.Time = TableData(:,1)';
if get(handles.checkboxElec,'Value')==1
    Start.Electricity = TableData(:,2)';
else Start.Electricity = 0*TableData(:,1)';
end
if get(handles.checkboxSteam,'Value')==1
    c = find(strcmp(ColNames,'Steam'));
    Start.Steam = TableData(:,c)';
else Start.Steam = TableData(:,1)';
end
if get(handles.checkboxHeat,'Value')==1
    c = find(strcmp(ColNames,'Heat'));
    Start.Heat = TableData(:,c)';
else Start.Heat = 0*TableData(:,1)';
end
if get(handles.checkboxCooling,'Value')==1
    c = find(strcmp(ColNames,'Cooling'));
    Start.Cooling = TableData(:,c)';
else Start.Cooling = 0*TableData(:,1)';
end
Start.Input = TableData(:,end)';
Generator.VariableStruct.Startup = Start;
plotStartup(handles)


% --- Executes on button press in pushbuttonShutDown.
function pushbuttonShutDown_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonShutDown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Generator EditTableHandle TableData
Finish = Generator.VariableStruct.Shutdown;
ColNames = {'Time'};
TableData =[];
TableData(:,1) = Finish.Time';
if get(handles.checkboxElec,'Value')==1
    ColNames(end+1) = {'Electricity'};
    TableData(:,end+1) = Finish.Electricity';
else Finish.Electricity = zeros(1,length(Finish.Time));
end
if get(handles.checkboxSteam,'Value')==1
    ColNames(end+1) = {'Steam'};
    TableData(:,end+1) = Finish.Steam';
else Finish.Steam = zeros(1,length(Finish.Time));
end
if get(handles.checkboxHeat,'Value')==1
    ColNames(end+1) = {'Heat'};
    TableData(:,end+1) = Finish.Heat';
else Finish.Heat = zeros(1,length(Finish.Time));
end
if get(handles.checkboxCooling,'Value')==1
    ColNames(end+1) = {'Cooling'};
    TableData(:,end+1) = Finish.Cooling;
else Finish.Cooling = zeros(1,length(Finish.Time));
end
TableData(:,end+1) = Finish.Input';
Text1 = 'Values represent normalized output 1 = 100%'; 
EditTable(ColNames,Text1);
waitfor(EditTableHandle)
Data = TableData;
Finish.Time = Data(:,1)';
if get(handles.checkboxElec,'Value')==1
    Finish.Electricity = Data(:,2)';
else Finish.Electricity = 0*Data(:,1)';
end
if get(handles.checkboxSteam,'Value')==1
    c = find(strcmp(ColNames,'Steam'));
    Finish.Steam = Data(:,c)';
else Finish.Steam = Data(:,1)';
end
if get(handles.checkboxHeat,'Value')==1
    c = find(strcmp(ColNames,'Heat'));
    Finish.Heat = Data(:,c)';
else Finish.Heat = 0*Data(:,1)';
end
if get(handles.checkboxCooling,'Value')==1
    c = find(strcmp(ColNames,'Cooling'));
    Finish.Cooling = Data(:,c)';
else Finish.Cooling = 0*Data(:,1)';
end
Finish.Input = Data(:,end)';
Generator.VariableStruct.Shutdown = Finish;
plotStartup(handles)

function plotStartup(handles)
global Generator
Start = Generator.VariableStruct.Startup;
Finish = Generator.VariableStruct.Shutdown;
axes(handles.axes3)
hold off
if get(handles.checkboxElec,'Value')==1
    plot(Start.Time,Start.Electricity,'k-')
    hold on
    plot(Finish.Time,Finish.Electricity,'k--')
end
if get(handles.checkboxSteam,'Value')==1
    plot(Start.Time,Start.Steam,'m-')
    hold on
    plot(Finish.Time,Finish.Steam,'m--')
end
if get(handles.checkboxHeat,'Value')==1
    plot(Start.Time,Start.Heat,'r-')
    hold on
    plot(Finish.Time,Finish.Heat,'r--')
end
if get(handles.checkboxCooling,'Value')==1
    plot(Start.Time,Start.Cooling,'b-')
    hold on
    plot(Finish.Time,Finish.Cooling,'b--')
end


% --- Executes on button press in pushbuttonSaveOnlyToPlant.
function pushbuttonSaveOnlyToPlant_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonSaveOnlyToPlant (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Plant SYSINDEX Generator

if nnz(nnz(Generator.Output.Electricity))~=get(handles.checkboxElec,'Value') || nnz(nnz(Generator.Output.Steam))~=get(handles.checkboxSteam,'Value') || nnz(nnz(Generator.Output.Heat))~=get(handles.checkboxHeat,'Value') || nnz(nnz(Generator.Output.Cooling))~=get(handles.checkboxCooling,'Value')
    errordlg('Output parameters do not match. Please edit output performance table.', 'Error Dialog');
    pushbuttonEditTable_Callback(hObject, eventdata, handles)
else
    if max(Generator.Output.Electricity)>0 && max(Generator.Output.Heat)>0
        Generator.Type = 'CHP Generator';
    elseif max(Generator.Output.Electricity)>0 
        Generator.Type = 'Electric Generator';
    elseif max(Generator.Output.Heat)>0
        Generator.Type = 'Heater';
    elseif max(Generator.Output.Cooling)>0 
        Generator.Type = 'Chiller';
        Generator.Size = 3.517.*get(handles.editSize,'Value');%conversion factor: cooling in kW = 3.517*cooling in tons
    elseif max(Generator.Output.Steam)>0 
        Generator.Type = 'Boiler';
    end
    Plant.Generator(SYSINDEX) = Generator;
    clear global Generator
    close(gcf)
end
% --- Executes on button press in pushbuttonCancel.
function pushbuttonCancel_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonCancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close(gcf)

% --- Executes on button press in pushbuttonSaveAs.
function pushbuttonSaveAs_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonSaveAs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global  Model_dir Generator SIMOhandle
Products = get(handles.checkboxElec,'Value')+get(handles.checkboxCooling,'Value')+get(handles.checkboxSteam,'Value')+get(handles.checkboxHeat,'Value');
r = size(Generator.VariableStruct.StateSpace.C);
if nnz(nnz(Generator.Output.Electricity))~=get(handles.checkboxElec,'Value') || nnz(nnz(Generator.Output.Steam))~=get(handles.checkboxSteam,'Value') || nnz(nnz(Generator.Output.Heat))~=get(handles.checkboxHeat,'Value') || nnz(nnz(Generator.Output.Cooling))~=get(handles.checkboxCooling,'Value')
    errordlg('Output parameters do not match. Please edit output performance table.', 'Error Dialog');
    pushbuttonEditTable_Callback(hObject, eventdata, handles)
elseif Products~=r(1)
    errordlg('Output matrix C dimensions do not match specified # of outputs. Please edit C.', 'Error Dialog');
    SIMOresponse;
    waitfor(SIMOhandle)
else
    if max(Generator.Output.Electricity)>0 && max(Generator.Output.Heat)>0
        savedir=fullfile(Model_dir,'component library','CHP Generator',Generator.Name);
        Generator.Type = 'CHP Generator';
    elseif max(Generator.Output.Electricity)>0 
        savedir=fullfile(Model_dir,'component library','Electric Generator',Generator.Name);
        Generator.Type = 'Electric Generator';
    elseif max(Generator.Output.Heat)>0
        savedir=fullfile(Model_dir,'component library','Heater',Generator.Name);
        Generator.Type = 'Heater';
    elseif max(Generator.Output.Cooling)>0 
        savedir=fullfile(Model_dir,'component library','Chiller',Generator.Name);
        Generator.Type = 'Chiller';
        Generator.Size = 3.517.*get(handles.editSize,'Value');%conversion factor: cooling in kW = 3.517*cooling in tons
    elseif max(Generator.Output.Steam)>0 
        savedir=fullfile(Model_dir,'component library','Boiler',Generator.Name);
        Generator.Type = 'Boiler';
    end
    if isempty(strfind(savedir,'.mat'))
        savedir = strcat(savedir,'.mat');
    end
    [f,p]=uiputfile(savedir,'Save As Generator Component');
    if f==0;return;end
    if isfield(Generator,'OpMatA')%if the model has already been run
        component = rmfield(Generator,'OpMatA');%remove OpMatA when saving the new component
    else
        component=Generator;    
    end
    save([p f],'component')
    pushbuttonSaveOnlyToPlant_Callback(hObject, eventdata, handles)
end

% --- Executes during object creation, after setting all properties.
function editDamp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editDamp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function editFreq_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editFreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbuttonComPorts.
function pushbuttonComPorts_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonComPorts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global EditCommHandle
CommunicationPorts();
waitfor(EditCommHandle)



function editRestartTime_Callback(hObject, eventdata, handles)
% hObject    handle to editRestartTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Generator
Generator.VariableStruct.RestartTime = str2double(get(hObject,'string'));

% --- Executes during object creation, after setting all properties.
function editRestartTime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editRestartTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
