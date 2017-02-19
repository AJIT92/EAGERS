function varargout = SetupElectricUtility(varargin)
% SETUPELECTRICUTILITY M-file for SetupElectricUtility.fig
%      SETUPELECTRICUTILITY, by itself, creates a new SETUPELECTRICUTILITY or raises the existing
%      singleton*.
%
%      H = SETUPELECTRICUTILITY returns the handle to a new SETUPELECTRICUTILITY or the handle to
%      the existing singleton*.
%
%      SETUPELECTRICUTILITY('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SETUPELECTRICUTILITY.M with the given input arguments.
%
%      SETUPELECTRICUTILITY('Property','Value',...) creates a new SETUPELECTRICUTILITY or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SetupElectricUtility_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SetupElectricUtility_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SetupElectricUtility

% Last Modified by GUIDE v2.5 18-May-2016 15:28:00

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SetupElectricUtility_OpeningFcn, ...
                   'gui_OutputFcn',  @SetupElectricUtility_OutputFcn, ...
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


% --- Executes just before SetupElectricUtility is made visible.
function SetupElectricUtility_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SetupElectricUtility (see VARARGIN)

% Choose default command line output for SetupElectricUtility
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
global Utility figHandle Timestamp
Utility = varargin{1};
Timestamp = varargin{2};
figHandle = handles.figure1;
set(handles.editName,'string',Utility.Name)
set(handles.uitable1,'data',Utility.SumRateTable)
set(handles.uitable2,'data',Utility.WinRateTable)
set(handles.popupmenuSummerMonth,'value',Utility.SumStartMonth)
m_d=[1 31; 2 28; 3 31; 4 30; 5 31; 6 30; 7 31; 8 31; 9 30; 10 31; 11 30; 12 31];
days=(1:m_d(Utility.SumStartMonth,2))';
set(handles.popupmenuSummerDay,'string',days,'value',Utility.SumStartDay)
set(handles.popupmenuWinterMonth,'value',Utility.WinStartMonth)
days=(1:m_d(Utility.WinStartMonth,2))';
set(handles.popupmenuWinterDay,'string',days,'value',Utility.WinStartDay)
set(handles.editSumEn1,'string',Utility.SumRates(1,1))
set(handles.editSumEn2,'string',Utility.SumRates(2,1))
set(handles.editSumEn3,'string',Utility.SumRates(3,1))
set(handles.editWinEn1,'string',Utility.WinRates(1,1))
set(handles.editWinEn2,'string',Utility.WinRates(2,1))
set(handles.editWinEn3,'string',Utility.WinRates(3,1))
set(handles.editSumDem1,'string',Utility.SumRates(1,2))
set(handles.editSumDem2,'string',Utility.SumRates(2,2))
set(handles.editSumDem3,'string',Utility.SumRates(3,2))
set(handles.editWinDem1,'string',Utility.WinRates(1,2))
set(handles.editWinDem2,'string',Utility.WinRates(2,2))
set(handles.editWinDem3,'string',Utility.WinRates(3,2))
set(handles.editMinThresh,'string',Utility.MinImportThresh)
%set(handles.editSellBack,'string',num2str(Utility.SellBackRate))
if ~isfield(Utility,'SellBackRate')
    Utility.SellBackRate = 0;
end
if ~isfield(Utility,'SellBackPerc')
    Utility.SellBackPerc =100;
end
set(handles.editSellbackPerc,'string',num2str(Utility.SellBackPerc))
if Utility.SellBackRate == 0
    set(handles.uipanelSellBack,'SelectedObject',handles.radiobuttonNoSellBack)
elseif Utility.SellBackRate == -1
    set(handles.uipanelSellBack,'SelectedObject',handles.radiobuttonReverseMeter)
else
    set(handles.uipanelSellBack,'SelectedObject',handles.radiobuttonPartialSellBack)
end


% --- Outputs from this function are returned to the command line.
function varargout = SetupElectricUtility_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


function editName_Callback(hObject, eventdata, handles)
% hObject    handle to editName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Utility
Utility.Name=get(hObject,'string');


% --- Executes during object creation, after setting all properties.
function editName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes when entered data in editable cell(s) in uitable1.
function uitable1_CellEditCallback(hObject, eventdata, handles)
global Utility
Utility.SumRateTable=get(handles.uitable1,'data');

% --- Executes when entered data in editable cell(s) in uitable2.
function uitable2_CellEditCallback(hObject, eventdata, handles)
global Utility
Utility.WinRateTable=get(handles.uitable2,'data');


% --- Executes on selection change in popupmenuSummerMonth.
function popupmenuSummerMonth_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuSummerMonth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Utility
Utility.SumStartMonth=get(hObject,'value');
m_d=[1 31; 2 28; 3 31; 4 30; 5 31; 6 30; 7 31; 8 31; 9 30; 10 31; 11 30; 12 31];
days=(1:m_d(Utility.SumStartMonth,2))';
if max(days)<Utility.SumStartDay
    set(handles.popupmenuSummerDay,'string',days,'value',max(days))
    popupmenuSummerDay_Callback(handles.popupmenuSummerDay,eventdata, handles)
else
    set(handles.popupmenuSummerDay,'string',Utility.SumStartDay)
end


% --- Executes during object creation, after setting all properties.
function popupmenuSummerMonth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuSummerMonth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in popupmenuSummerDay.
function popupmenuSummerDay_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuSummerDay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Utility
Utility.SumStartDay=get(hObject,'value');

% --- Executes during object creation, after setting all properties.
function popupmenuSummerDay_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuSummerDay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in popupmenuWinterMonth.
function popupmenuWinterMonth_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuWinterMonth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Utility
Utility.WinStartMonth=get(hObject,'value');
m_d=[1 31; 2 28; 3 31; 4 30; 5 31; 6 30; 7 31; 8 31; 9 30; 10 31; 11 30; 12 31];
days=(1:m_d(Utility.WinStartMonth,2))';
if max(days)<Utility.WinStartDay
    set(handles.popupmenuWinterDay,'string',days,'value',max(days))
    popupmenuWinterDay_Callback(handles.popupmenuWinterDay,eventdata, handles)
else
    set(handles.popupmenuWinterDay,'string',Utility.WinStartDay)
end

% --- Executes during object creation, after setting all properties.
function popupmenuWinterMonth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuWinterMonth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in popupmenuWinterDay.
function popupmenuWinterDay_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuWinterDay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Utility
Utility.WinStartDay=get(hObject,'value');

% --- Executes during object creation, after setting all properties.
function popupmenuWinterDay_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuWinterDay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function editWinEn1_Callback(hObject, eventdata, handles)
% hObject    handle to editWinEn1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Utility
Utility.WinRates(1,1)=str2double(get(hObject,'string'));

% --- Executes during object creation, after setting all properties.
function editWinEn1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editWinEn1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editWinEn2_Callback(hObject, eventdata, handles)
% hObject    handle to editWinEn2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Utility
Utility.WinRates(2,1)=str2double(get(hObject,'string'));

% --- Executes during object creation, after setting all properties.
function editWinEn2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editWinEn2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editWinEn3_Callback(hObject, eventdata, handles)
% hObject    handle to editWinEn3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Utility
Utility.WinRates(3,1)=str2double(get(hObject,'string'));

% --- Executes during object creation, after setting all properties.
function editWinEn3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editWinEn3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editSumEn1_Callback(hObject, eventdata, handles)
% hObject    handle to editSumEn1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Utility
Utility.SumRates(1,1)=str2double(get(hObject,'string'));


% --- Executes during object creation, after setting all properties.
function editSumEn1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editSumEn1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editSumEn2_Callback(hObject, eventdata, handles)
% hObject    handle to editSumEn2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Utility
Utility.SumRates(2,1)=str2double(get(hObject,'string'));

% --- Executes during object creation, after setting all properties.
function editSumEn2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editSumEn2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editSumEn3_Callback(hObject, eventdata, handles)
% hObject    handle to editSumEn3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Utility
Utility.SumRates(3,1)=str2double(get(hObject,'string'));

% --- Executes during object creation, after setting all properties.
function editSumEn3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editSumEn3 (see GCBO)
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
global Plant SYSINDEX Utility
if handles.radiobuttonReverseMeter.Value
    Utility.SellBackRate = -1;
elseif handles.radiobuttonPartialSellBack.Value
    Utility.SellBackPerc = str2double(get(handles.editSellBack,'String'))/100;
    Utility.SellBackRate = 0;
else Utility.SellBackRate = 0;
    Utility.SellBackPerc = 0;
end
Plant.Generator(SYSINDEX).VariableStruct = Utility;
close(gcf)

% --- Executes on button press in pushbuttonSaveAs.
function pushbuttonSaveAs_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonSaveAs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% x=get(handles.uitable1,'data');
global Utility Model_dir
Utility.MinImportThresh = str2double(get(handles.editMinThresh,'String'));
[f,p]=uiputfile(fullfile(Model_dir,'System Library','Utility',strcat(char(Utility.Name),'.mat')),'Save As Grid Component...');
if f==0;return;end
if isfield(Utility,'OpMatA')%if the model has already been run
    component = rmfield(Utility,'OpMatA');%remove OpMatA when saving the new component
else
    component=Utility;    
end
save([p f],'component')
pushbuttonSaveOnlyToProject_Callback(hObject, eventdata, handles)


function editWinDem1_Callback(hObject, eventdata, handles)
% hObject    handle to editWinDem1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Utility
Utility.WinRates(1,2)=str2double(get(hObject,'string'));


% --- Executes during object creation, after setting all properties.
function editWinDem1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editWinDem1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editWinDem2_Callback(hObject, eventdata, handles)
% hObject    handle to editWinDem2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Utility
Utility.WinRates(2,2)=str2double(get(hObject,'string'));

% --- Executes during object creation, after setting all properties.
function editWinDem2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editWinDem2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editWinDem3_Callback(hObject, eventdata, handles)
% hObject    handle to editWinDem3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Utility
Utility.WinRates(3,2)=str2double(get(hObject,'string'));


% --- Executes during object creation, after setting all properties.
function editWinDem3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editWinDem3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editSumDem1_Callback(hObject, eventdata, handles)
% hObject    handle to editSumDem1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Utility
Utility.SumRates(1,2)=str2double(get(hObject,'string'));

% --- Executes during object creation, after setting all properties.
function editSumDem1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editSumDem1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editSumDem2_Callback(hObject, eventdata, handles)
% hObject    handle to editSumDem2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Utility
Utility.SumRates(2,2)=str2double(get(hObject,'string'));


% --- Executes during object creation, after setting all properties.
function editSumDem2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editSumDem2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editSumDem3_Callback(hObject, eventdata, handles)
% hObject    handle to editSumDem3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Utility
Utility.SumRates(3,2)=str2double(get(hObject,'string'));


% --- Executes during object creation, after setting all properties.
function editSumDem3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editSumDem3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editMinThresh_Callback(hObject, eventdata, handles)
% hObject    handle to editMinThresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Utility
Utility.MinImportThresh=str2double(get(hObject,'string'));

% --- Executes during object creation, after setting all properties.
function editMinThresh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editMinThresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editSellbackPerc_Callback(hObject, eventdata, handles)
% hObject    handle to editSellbackPerc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Utility
Utility.SellBackRate =str2double(get(hObject,'string'))/100;

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
global Utility
switch get(eventdata.NewValue,'Tag')
    case 'radiobuttonNoSellBack'
        Utility.SellBackRate = 0;
        set(handles.editSellBack,'string','0.00')
    case 'radiobuttonPartialSellBack'
        Utility.SellBackRate = str2double(get(handles.editSellBack,'String'))/100;
        set(handles.editSellBack,'string',num2str(Utility.SellBackRate*100))
    case 'radiobuttonReverseMeter'
        Utility.SellBackRate = -1;
        set(handles.editSellBack,'string','-1')
end


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
delete(hObject);
