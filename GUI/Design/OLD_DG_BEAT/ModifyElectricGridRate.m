function varargout = ModifyElectricGridRate(varargin)
% MODIFYELECTRICGRIDRATE M-file for ModifyElectricGridRate.fig
%      MODIFYELECTRICGRIDRATE, by itself, creates a new MODIFYELECTRICGRIDRATE or raises the existing
%      singleton*.
%
%      H = MODIFYELECTRICGRIDRATE returns the handle to a new MODIFYELECTRICGRIDRATE or the handle to
%      the existing singleton*.
%
%      MODIFYELECTRICGRIDRATE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MODIFYELECTRICGRIDRATE.M with the given input arguments.
%
%      MODIFYELECTRICGRIDRATE('Property','Value',...) creates a new MODIFYELECTRICGRIDRATE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ModifyElectricGridRate_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ModifyElectricGridRate_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ModifyElectricGridRate

% Last Modified by GUIDE v2.5 07-Apr-2014 11:23:58

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ModifyElectricGridRate_OpeningFcn, ...
                   'gui_OutputFcn',  @ModifyElectricGridRate_OutputFcn, ...
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


% --- Executes just before ModifyElectricGridRate is made visible.
function ModifyElectricGridRate_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ModifyElectricGridRate (see VARARGIN)

% Choose default command line output for ModifyElectricGridRate
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

global Project COMPONENT Model_dir
COMPONENT=Project.Utilities.Grid;%Rate structure

compdir=fullfile(Model_dir,'System Library', 'Grid');
files=dir(fullfile(compdir,'*.mat'));
list=strrep({files.name},'.mat','');
val = find(strcmp(COMPONENT.Name,list),1);
if isempty(val)
    val = 1;
end
set(handles.popupmenuUtility,'string',list,'value',val)
listUser = {'Residential'; 'Commercial';'Industrial';};
set(handles.popupmenuUser,'string',listUser,'value',2)
stateName = {'Alabama';'Alaska';'Arizona';'Arkansas';'California';'Colorado';'Connecticut';'Delaware';'Florida';'Georgia';
             'Hawaii';'Idaho';'Illinois';'Indiana';'Iowa';'Kansas';'Kentucky';'Louisiana';'Maine';'Maryland';
             'Massachusetts';'Michigan';'Minnesota';'Mississippi';'Missouri';'Montana';'Nebraska';'Nevada';'NewHampshire';'NewJersey';
             'NewMexico';'NewYork';'NorthCarolina';'NorthDakota';'Ohio';'Oklahoma';'Oregon';'Pennsylvania';'RhodeIsland';'SouthCarolina';
             'SouthDakota';'Tennessee';'Texas';'Utah';'Vermont';'Virginia';'Washington';'WestVirginia';'Wisconsin';'Wyoming';};
val = find(strcmp(Project.State,stateName),1);
if isempty(val)
    val = 1;
end
set(handles.popupmenuState,'string',stateName,'value',val)
updateGUI(handles)   
plotElecRate(handles)


function updateGUI(handles)
global COMPONENT
set(handles.editName,'string',COMPONENT.Name)
set(handles.editDescription,'string',COMPONENT.Description)

set(handles.uitable1,'data',COMPONENT.summerRateTable)
set(handles.uitable2,'data',COMPONENT.winterRateTable)
datademand=COMPONENT.useRatesDescription;
for i=1:length(COMPONENT.demandCharges)
    datademand{i,2}=COMPONENT.demandCharges(i);
end
set(handles.uitableDemand,'data',datademand)
set(handles.popupmenuSummerFromMonth,'value',COMPONENT.summerStartDate(1))
m_d=[1 31; 2 28; 3 31; 4 30; 5 31; 6 30; 7 31; 8 31; 9 30; 10 31; 11 30; 12 31];
days=(1:m_d(COMPONENT.summerStartDate(1),2))';
set(handles.popupmenuSummerFromDay,'string',days,'value',COMPONENT.summerStartDate(2))
set(handles.popupmenuSummerToMonth,'value',COMPONENT.summerEndDate(1))
days=(1:m_d(COMPONENT.summerEndDate(1),2))';
set(handles.popupmenuSummerToDay,'string',days,'value',COMPONENT.summerEndDate(2))
clear days
set(handles.editSummer1,'string',COMPONENT.summerPrice(1))
set(handles.editSummer2,'string',COMPONENT.summerPrice(2))
set(handles.editSummer3,'string',COMPONENT.summerPrice(3))
set(handles.editWinter1,'string',COMPONENT.winterPrice(1))
set(handles.editWinter2,'string',COMPONENT.winterPrice(2))
set(handles.editWinter3,'string',COMPONENT.winterPrice(3))
set(handles.editSellBack,'string',num2str(COMPONENT.SellBackRate))
if ~isfield(COMPONENT,'SellBackPerc')
    COMPONENT.SellBackPerc =100;
end
set(handles.editSellbackPerc,'string',num2str(COMPONENT.SellBackPerc))
if COMPONENT.SellBackRate == 0
    set(handles.uipanelSellBack,'SelectedObject',handles.radiobuttonNoSellBack)
elseif COMPONENT.SellBackRate == -1
    set(handles.uipanelSellBack,'SelectedObject',handles.radiobuttonReverseMeter)
else
    set(handles.uipanelSellBack,'SelectedObject',handles.radiobuttonFixedSellBack)
end


% --- Outputs from this function are returned to the command line.
function varargout = ModifyElectricGridRate_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes when entered data in editable cell(s) in uitable1.
function uitable1_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to uitable1 (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
global COMPONENT
COMPONENT.summerRateTable=get(handles.uitable1,'data');


% --- Executes on button press in pushbuttonSaveAs.
function pushbuttonSaveAs_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonSaveAs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Project COMPONENT Model_dir
[f,p]=uiputfile('*.mat','Save As Electric Rate Structure',fullfile(Model_dir,'System Library','Grid'));
if f==0;return;end
component=COMPONENT;
save([p f],'component')
Project.Utilities.Grid=COMPONENT;
clear global COMPONENT 
close(gcf)
uiresume

% --- Executes on selection change in popupmenuSummerFromMonth.
function popupmenuSummerFromMonth_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuSummerFromMonth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenuSummerFromMonth contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenuSummerFromMonth
global COMPONENT
COMPONENT.summerStartDate(1)=get(hObject,'value');
m_d=[1 31; 2 28; 3 31; 4 30; 5 31; 6 30; 7 31; 8 31; 9 30; 10 31; 11 30; 12 31];
days=(1:m_d(COMPONENT.summerStartDate(1),2))';
if max(days)<COMPONENT.summerStartDate(2)
    set(handles.popupmenuSummerFromDay,'string',days,'value',max(days))
    popupmenuSummerFromDay_Callback(handles.popupmenuSummerFromDay,eventdata, handles)
else
    set(handles.popupmenuSummerFromDay,'string',days)
end


% --- Executes during object creation, after setting all properties.
function popupmenuSummerFromMonth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuSummerFromMonth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on selection change in popupmenuSummerFromDay.
function popupmenuSummerFromDay_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuSummerFromDay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global COMPONENT
COMPONENT.summerStartDate(2)=get(hObject,'value');


% --- Executes during object creation, after setting all properties.
function popupmenuSummerFromDay_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuSummerFromDay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenuSummerToMonth.
function popupmenuSummerToMonth_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuSummerToMonth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global COMPONENT
COMPONENT.summerEndDate(1)=get(hObject,'value');
m_d=[1 31; 2 28; 3 31; 4 30; 5 31; 6 30; 7 31; 8 31; 9 30; 10 31; 11 30; 12 31];
days=(1:m_d(COMPONENT.summerEndDate(1),2))';
if max(days)<COMPONENT.summerEndDate(2)
    set(handles.popupmenuSummerToDay,'string',days,'value',max(days))
    popupmenuSummerToDay_Callback(handles.popupmenuSummerToDay,eventdata, handles)
else
    set(handles.popupmenuSummerToDay,'string',days)
end

% --- Executes during object creation, after setting all properties.
function popupmenuSummerToMonth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuSummerToMonth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenuSummerToDay.
function popupmenuSummerToDay_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuSummerToDay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global COMPONENT
COMPONENT.summerEndDate(2)=get(hObject,'value');

% --- Executes during object creation, after setting all properties.
function popupmenuSummerToDay_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuSummerToDay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function editWinter1_Callback(hObject, eventdata, handles)
% hObject    handle to editWinter1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global COMPONENT
COMPONENT.winterPrice(1)=str2double(get(hObject,'string'));



% --- Executes during object creation, after setting all properties.
function editWinter1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editWinter1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editWinter2_Callback(hObject, eventdata, handles)
% hObject    handle to editWinter2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global COMPONENT
COMPONENT.winterPrice(2)=str2double(get(hObject,'string'));

% --- Executes during object creation, after setting all properties.
function editWinter2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editWinter2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editWinter3_Callback(hObject, eventdata, handles)
% hObject    handle to editWinter3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global COMPONENT
COMPONENT.winterPrice(3)=str2double(get(hObject,'string'));

% --- Executes during object creation, after setting all properties.
function editWinter3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editWinter3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editSummer1_Callback(hObject, eventdata, handles)
% hObject    handle to editSummer1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global COMPONENT
COMPONENT.summerPrice(1)=str2double(get(hObject,'string'));


% --- Executes during object creation, after setting all properties.
function editSummer1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editSummer1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editSummer2_Callback(hObject, eventdata, handles)
% hObject    handle to editSummer2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global COMPONENT
COMPONENT.summerPrice(2)=str2double(get(hObject,'string'));

% --- Executes during object creation, after setting all properties.
function editSummer2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editSummer2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editSummer3_Callback(hObject, eventdata, handles)
% hObject    handle to editSummer3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global COMPONENT
COMPONENT.summerPrice(3)=str2double(get(hObject,'string'));

% --- Executes during object creation, after setting all properties.
function editSummer3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editSummer3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when entered data in editable cell(s) in uitable2.
function uitable2_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to uitable2 (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
global COMPONENT
COMPONENT.winterRateTable=get(handles.uitable2,'data');


% --- Executes on button press in pushbuttonCancel.
function pushbuttonCancel_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonCancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clear global COMPONENT
close(gcf)


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
clear global COMPONENT
delete(hObject);

function editDescription_Callback(hObject, eventdata, handles)
% hObject    handle to editDescription (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global COMPONENT
COMPONENT.Description=get(hObject,'string');

% --- Executes during object creation, after setting all properties.
function editDescription_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editDescription (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editName_Callback(hObject, eventdata, handles)
% hObject    handle to editName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global COMPONENT
COMPONENT.Name=get(hObject,'string');

% --- Executes during object creation, after setting all properties.
function editName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes when entered data in editable cell(s) in uitableDemand.
function uitableDemand_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to uitableDemand (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
global COMPONENT
data=get(hObject,'data');
COMPONENT.demandCharges=[data{:,2}]';


% --- Executes on button press in pushbuttonSaveOnlyToProject.
function pushbuttonSaveOnlyToProject_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonSaveOnlyToProject (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Project COMPONENT 
Project.Utilities.Grid=COMPONENT;
clear global COMPONENT
close(gcf)
uiresume

function editSellBack_Callback(hObject, eventdata, handles)
% hObject    handle to editSellBack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global COMPONENT 
COMPONENT.SellBackRate = str2double(get(hObject,'String'));

% --- Executes during object creation, after setting all properties.
function editSellBack_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editSellBack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editSellbackPerc_Callback(hObject, eventdata, handles)
% hObject    handle to editSellbackPerc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global COMPONENT
COMPONENT.SellBackPerc =str2double(get(hObject,'string'));

% --- Executes during object creation, after setting all properties.
function editSellbackPerc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editSellbackPerc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes when selected object is changed in uipanelChillType.
function uipanelSellBack_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanelChillType 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
global COMPONENT
switch get(eventdata.NewValue,'Tag')
    case 'radiobuttonNoSellBack'
        COMPONENT.SellBackRate = 0;
        set(handles.editSellBack,'string','0.00')
    case 'radiobuttonFixedSellBack'
        if COMPONENT.SellBackRate < 0;
            COMPONENT.SellBackRate = 0;
        else
            COMPONENT.SellBackRate = str2double(get(handles.editSellBack,'String'));
        end
        set(handles.editSellBack,'string',num2str(COMPONENT.SellBackRate))
    case 'radiobuttonReverseMeter'
        COMPONENT.SellBackRate = -1;
        set(handles.editSellBack,'string','-1')
end


% --- Executes on selection change in popupmenuUtility.
function popupmenuUtility_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuUtility (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global COMPONENT Model_dir
list = get(handles.popupmenuUtility,'string');
val = get(handles.popupmenuUtility,'value');
compdir=fullfile(Model_dir,'System Library', 'Grid');
load(strcat(compdir,filesep,list{val})) %load component structure
COMPONENT = component;
set(handles.popupmenuUtility,'string',list,'value',val)
updateGUI(handles)


% --- Executes during object creation, after setting all properties.
function popupmenuUtility_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuUtility (see GCBO)
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
global Project
list = get(handles.popupmenuState,'string');
val = get(handles.popupmenuState,'value');
Project.State = list(val);
plotElecRate(handles)

% --- Executes during object creation, after setting all properties.
function popupmenuState_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuState (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function plotElecRate(handles)
global Project Model_dir
str = get(handles.popupmenuUser,'string');
val = get(handles.popupmenuUser,'Value');
user = str{val};
StateNames = get(handles.popupmenuState,'string');
stateAbrev = {'AL';'AK';'AZ';'AR';'CA';'CO';'CT';'DE';'FL';'GA';'HI';'ID';'IL';'IN';'IA';'KS';'KY';'LA';'ME';'MD';'MA';'MI';'MN';'MS';'MO';
              'MT';'NE';'NV';'NH';'NJ';'NM';'NY';'NC';'ND';'OH';'OK';'OR';'PA';'RI';'SC';'SD';'TN';'TX';'UT';'VT';'VA';'WA';'WV';'WI';'WY';};
I = find(strcmp(Project.State,StateNames),1);
state = char(stateAbrev(I));
Projection = ProjectUtilityCost('electric',user,10,state);
load(fullfile(Model_dir,'System Library', 'Grid','RateData','ElecRate'))
Date = ElecRate.Date;
Residential = ElecRate.(state).Residential;
Commercial = ElecRate.(state).Commercial;
Industrial = ElecRate.(state).Industrial;
axes(handles.axes1)
cla reset
set(gca,'tag','axes1')
% plot(Date,AllSector,'r')
hold on
plot(Date,Residential,'g')
plot(Date,Commercial,'b')
plot(Date,Industrial,'m')

for j = 1:1:120
    Date2(j) =  datenum(2014, j,01);
end
plot(Date2,Projection,'k')
DatePlot = [Date(1:24:end) Date2(1:24:end)];
set(gca,'xtick',DatePlot)
datetick('x','yy','keepticks')
xlabel('Year')
ylabel('Electricity Cost [cents/kW]')
legend('Residential' ,'Commercial','Industrial','Projecion', 'Location','NorthWest')
title(char(Project.State))


% --- Executes on selection change in popupmenuUser.
function popupmenuUser_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuUser (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
plotElecRate(handles);

% --- Executes during object creation, after setting all properties.
function popupmenuUser_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuUser (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

