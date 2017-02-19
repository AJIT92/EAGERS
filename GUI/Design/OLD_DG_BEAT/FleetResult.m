function varargout = FleetResult(varargin)
% FLEETRESULT M-file for FleetResult.fig
%      FLEETRESULT, by itself, creates a new FLEETRESULT or raises the existing
%      singleton*.
%
%      H = FLEETRESULT returns the handle to a new FLEETRESULT or the handle to
%      the existing singleton*.
%
%      FLEETRESULT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FLEETRESULT.M with the given input arguments.
%
%      FLEETRESULT('Property','Value',...) creates a new FLEETRESULT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before FleetResult_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to FleetResult_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help FleetResult

% Last Modified by GUIDE v2.5 22-Jan-2014 12:14:42

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @FleetResult_OpeningFcn, ...
                   'gui_OutputFcn',  @FleetResult_OutputFcn, ...
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


% --- Executes just before FleetResult is made visible.
function FleetResult_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to FleetResult (see VARARGIN)

% Choose default command line output for FleetResult
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

global Model_dir 
resultsdir=fullfile(Model_dir, 'results'); 
x=dir(fullfile(resultsdir,'*.mat'));
listTemp=strrep({x.name},'.mat','');
i = 0;
for j = 1:1:length(listTemp')
    if  strncmp(listTemp(j),'Fleet',5)
        i = i+1;
        list(i) = strrep(listTemp(j),'Fleet','');
    end
end
if ~isempty(find(strcmp('Project1',list)));
    currentProject = find(strcmp('Project1',list));
else currentProject=1;
end
set(handles.popupmenuResultsFile,'string',list,'value',currentProject)
list2 = {'NPV'; 'CO2'; 'NOx'; 'SO2';};
set(handles.popupmenuNPVemission,'string',list2,'value',1)
popupmenuResultsFile_Callback(handles.popupmenuResultsFile, eventdata, handles)


% --- Outputs from this function are returned to the command line.
function varargout = FleetResult_OutputFcn(hObject, eventdata, handles) 
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
global Model_dir Fleet 
str=get(hObject,'string');
val=get(hObject,'value'); 
resultsdir=fullfile(Model_dir, 'results');
load(fullfile(resultsdir, strcat('Fleet',str{val}))) % loads Fleet
updateTable(handles)

%% fill in table
function updateTable(handles)
global TableData Fleet
n = length(Fleet.BuildName);
TableData ={};
TableData(:,1) = Fleet.BuildName;
for i = 1:1:n
    TableData(i,2) = cellstr(num2str(Fleet.InstalledCHPcost(i)));
    NPVsave = sum(Fleet.NPVbaseline(1:5,i)) - sum(Fleet.NPVdispatch(1:5,i));
    NPVpercSave = NPVsave/sum(Fleet.NPVbaseline(1:5,i))*100;
    TableData(i,3) = cellstr(num2str(NPVsave));
    TableData(i,4) = cellstr(num2str(NPVpercSave));
    GHGsave = sum(Fleet.Summary.CO2baseline(i,:))-sum(Fleet.Summary.CO2dispatch(i,:));
    TableData(i,5) = cellstr(num2str(GHGsave));
end

set(handles.uitableResult,'Data',TableData)

updatePlots(handles)

%% Plot
function updatePlots(handles)
global Fleet BuildSelect
val = BuildSelect;
if isempty(val)
    val =1;
end
val2 = get(handles.popupmenuNPVemission,'value');
str=get(handles.popupmenuNPVemission,'string');
currentstr=str{val2};

axes(handles.axes1)
cla reset
set(gca,'tag','axes1') 

switch currentstr
     case 'NPV'
        bar([Fleet.NPVbaseline(1:5,val)';Fleet.NPVdispatch(1:5,val)';],'stacked')
        ylabel('Fuel, Grid, and O&M/Finance Costs ($)')
        set(gca,'XTickLabel',{'Baseline','System'})
        xlim([0.5 2.5])
        legend('Demand', 'Grid' ,'Fuel', 'O&M','Financing','Location','NorthEastOutside')
	case 'CO2'
        bar([Fleet.Summary.CO2baseline(val,1:3);Fleet.Summary.CO2dispatch(val,1:3);],'stacked')
        ylabel('Emissions of CO2 (lbs)')
        set(gca,'XTickLabel',{'Baseline','System'})
        xlim([0.5 2.5])
        legend('Grid', 'CHP' ,'Boiler', 'Location','NorthEastOutside')
    case 'NOx'
        bar([Fleet.Summary.NOxbaseline(val,1:3);Fleet.Summary.NOxdispatch(val,1:3);],'stacked')
        ylabel('Emissions of  NOx (lbs)')
        set(gca,'XTickLabel',{'Baseline','System'})
        xlim([0.5 2.5])
        legend('Grid', 'CHP' ,'Boiler', 'Location','NorthEastOutside')
    case 'SO2'
        bar([Fleet.Summary.SO2baseline(val,1:3);Fleet.Summary.SO2dispatch(val,1:3);],'stacked')
        ylabel('Emissions of SO2 (lbs)')
        set(gca,'XTickLabel',{'Baseline','System'})
        xlim([0.5 2.5])
        legend('Grid', 'CHP' ,'Boiler', 'Location','NorthEastOutside')               
end


% --- Executes during object creation, after setting all properties.
function popupmenuResultsFile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuResultsFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pushbuttonSortCost.
function pushbuttonSortCost_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonSortCost (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global TableData
A = str2double(TableData(:,2));
if ~issorted(A)
    [Y, I] = sort(A,1,'descend');
else
    [Y, I] = sort(A,1,'ascend');
end

for i = 1:1:5
    Temporary(:,1) = TableData(I,i);
end
TableData = Temporary;
set(handles.uitableResult,'Data',TableData)

% --- Executes on button press in pushbuttonSortNPV.
function pushbuttonSortNPV_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonSortNPV (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global TableData
A = str2double(TableData(:,3));
if ~issorted(A)
    [Y, I] = sort(A,1,'descend');
else
    [Y, I] = sort(A,1,'ascend');
end

for i = 1:1:5
    Temporary(:,1) = TableData(I,i);
end
TableData = Temporary;
set(handles.uitableResult,'Data',TableData)

% --- Executes on button press in pushbuttonSortPercNPV.
function pushbuttonSortPercNPV_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonSortPercNPV (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global TableData
A = str2double(TableData(:,4));
if ~issorted(A)
    [Y, I] = sort(A,1,'descend');
else
    [Y, I] = sort(A,1,'ascend');
end

for i = 1:1:5
    Temporary(:,1) = TableData(I,i);
end
TableData = Temporary;
set(handles.uitableResult,'Data',TableData)

% --- Executes on button press in pushbuttonSortGHG.
function pushbuttonSortGHG_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonSortGHG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global TableData
A = str2double(TableData(:,5));
if ~issorted(A)
    [Y, I] = sort(A,1,'descend');
else
    [Y, I] = sort(A,1,'ascend');
end

for i = 1:1:5
    Temporary(:,1) = TableData(I,i);
end
TableData = Temporary;
set(handles.uitableResult,'Data',TableData)

% --- Executes on selection change in popupmenuNPVemission.
function popupmenuNPVemission_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuNPVemission (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
updatePlots(handles)

% --- Executes during object creation, after setting all properties.
function popupmenuNPVemission_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuNPVemission (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes when selected cell(s) is changed in uitableResult.
function uitableResult_CellSelectionCallback(hObject, eventdata, handles)
% hObject    handle to uitableResult (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)
global BuildSelect
BuildSelect = eventdata.Indices(1);
updatePlots(handles)


