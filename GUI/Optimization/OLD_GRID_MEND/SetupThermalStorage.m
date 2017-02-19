function varargout = SetupThermalStorage(varargin)
% SETUPTHERMALSTORAGE MATLAB code for SetupThermalStorage.fig
%      SETUPTHERMALSTORAGE, by itself, creates a new SETUPTHERMALSTORAGE or raises the existing
%      singleton*.
%
%      H = SETUPTHERMALSTORAGE returns the handle to a new SETUPTHERMALSTORAGE or the handle to
%      the existing singleton*.
%
%      SETUPTHERMALSTORAGE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SETUPTHERMALSTORAGE.M with the given input arguments.
%
%      SETUPTHERMALSTORAGE('Property','Value',...) creates a new SETUPTHERMALSTORAGE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SetupThermalStorage_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SetupThermalStorage_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SetupThermalStorage

% Last Modified by GUIDE v2.5 14-May-2015 10:18:27

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SetupThermalStorage_OpeningFcn, ...
                   'gui_OutputFcn',  @SetupThermalStorage_OutputFcn, ...
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


% --- Executes just before SetupThermalStorage is made visible.
function SetupThermalStorage_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SetupThermalStorage (see VARARGIN)

% Choose default command line output for SetupThermalStorage
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

global Plant SYSINDEX figHandle Generator
figHandle = handles.figure1;
Generator = Plant.Generator(SYSINDEX);
TES=Plant.Generator(SYSINDEX).VariableStruct;

set(handles.editName,'string',Generator.Name)
% set(handles.textType,'string',COMPONENT.Type)
set(handles.editSize,'string',Generator.Size)
if strcmp(TES.EnStoreType, 'ColdTES')
    set(handles.uipanelEnStoreType,'SelectedObject', handles.radiobuttonColdTES)
elseif strcmp(TES.EnStoreType, 'HotTES')
    set(handles.uipanelEnStoreType,'SelectedObject', handles.radiobuttonHotTES)
end

%TES specific
set(handles.editSizeLiter,'string',TES.SizeLiter)
set(handles.editTcold,'string',TES.Tcold)
set(handles.editThot,'string',TES.Thot)
set(handles.editLoss,'string',TES.SelfDischarge*100*24)
set(handles.editFillRate,'string',TES.FillRate)
set(handles.editFillRatePerc,'string',TES.FillRatePerc)
set(handles.editDischRate,'string',TES.DischRate)
set(handles.editDischRatePerc,'string',TES.DischRatePerc)
set(handles.editDischargeEff,'string',TES.DischargeEff*100)
set(handles.editChargeEff,'string',TES.ChargeEff*100)

fillRatesStatic(handles)


% --- Outputs from this function are returned to the command line.
function varargout = SetupThermalStorage_OutputFcn(hObject, eventdata, handles) 
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
global Generator
kWh2Liter = 3600/4.186; %3600 seconds / 4.186 kJ/kg*K 
Generator.Size=str2double(get(hObject,'string'));
Generator.VariableStruct.SizeLiter = Generator.Size*kWh2Liter/(Generator.VariableStruct.Thot-Generator.VariableStruct.Tcold);
set(handles.editSizeLiter,'string',Generator.VariableStruct.SizeLiter)

% --- Executes during object creation, after setting all properties.
function editSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function editName_Callback(hObject, eventdata, handles)
% hObject    handle to editName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Generator
Generator.Name=get(hObject,'string');


% --- Executes during object creation, after setting all properties.
function editName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbuttonSaveOnlyToPlant.
function pushbuttonSaveOnlyToPlant_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonSaveOnlyToPlant (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Plant SYSINDEX Generator
Plant.Generator(SYSINDEX) = Generator;
clear global Generator
close(gcf)

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
global Generator Model_dir
if isfield(Generator,'OpMatA')%if the model has already been run
    component = rmfield(Generator,'OpMatA');%remove OpMatA when saving the new component
else
    component=Generator;    
end
savedir=fullfile(Model_dir,'System Library','Thermal Storage',strcat(Generator.Name,'.mat'));
[f,p]=uiputfile(savedir,'Save As Energy Storage Component');
if f==0;return;end
save([p f],'component')
pushbuttonSaveOnlyToPlant_Callback(hObject, eventdata, handles)



function editSizeLiter_Callback(hObject, eventdata, handles)
% hObject    handle to editSizeLiter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Generator
Generator.VariableStruct.SizeLiter=str2double(get(hObject,'string'));
kWh2Liter = 3600/4.186; %3600 seconds / 4.186 kJ/kg*K
Generator.Size = Generator.VariableStruct.SizeLiter/kWh2Liter*(Generator.VariableStruct.Thot-Generator.VariableStruct.Tcold);
set(handles.editSize,'string',Generator.Size)


% --- Executes during object creation, after setting all properties.
function editSizeLiter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editSizeLiter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editTcold_Callback(hObject, eventdata, handles)
% hObject    handle to editTcold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Generator
Generator.VariableStruct.Tcold=str2double(get(hObject,'string'));
kWh2Liter = 3600/4.186; %3600 seconds / 4.186 kJ/kg*K 
Generator.Size = Generator.VariableStruct.SizeLiter/kWh2Liter*(Generator.VariableStruct.Thot-Generator.VariableStruct.Tcold);
set(handles.editSize,'string',Generator.Size)


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
global Generator
Generator.VariableStruct.Thot=str2double(get(hObject,'string'));
kWh2Liter = 3600/4.186; %3600 seconds / 4.186 kJ/kg*K 
Generator.Size = Generator.VariableStruct.SizeLiter/kWh2Liter*(Generator.VariableStruct.Thot-Generator.VariableStruct.Tcold);
set(handles.editSize,'string',Generator.Size)


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
global Generator
Generator.VariableStruct.FillRate = str2double(get(hObject,'string'));% Fill rate in liter/min
Generator.VariableStruct.FillRatePerc = Generator.VariableStruct.FillRate*60/Generator.VariableStruct.SizeLiter*100; % Fill rate in (%/hr)
set(handles.editFillRatePerc,'string',Generator.VariableStruct.FillRatePerc)


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
global Generator
Generator.VariableStruct.FillRatePerc = str2double(get(hObject,'string'));
Generator.VariableStruct.FillRate = Generator.VariableStruct.SizeLiter*Generator.VariableStruct.FillRatePerc/100/60; % Fill rate in liter/min
set(handles.editFillRate,'string',Generator.VariableStruct.FillRate)

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
global Generator
Generator.VariableStruct.DischRate = str2double(get(hObject,'string'));% disch rate in liter/min
Generator.VariableStruct.DischRatePerc = Generator.VariableStruct.DischRate*60/Generator.VariableStruct.SizeLiter*100; % disch rate in (%/hr)
set(handles.editDischRatePerc,'string',Generator.VariableStruct.DischRatePerc)

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
global Generator
Generator.VariableStruct.DischRatePerc = str2double(get(hObject,'string'));% disch rate in (%/hr)
Generator.VariableStruct.DischRate = Generator.VariableStruct.SizeLiter*Generator.VariableStruct.DischRatePerc/100/60; % disch rate in liter/min
set(handles.editDischRate,'string',Generator.VariableStruct.DischRate)

% --- Executes during object creation, after setting all properties.
function editDischRatePerc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editDischRatePerc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editChargeEff_Callback(hObject, eventdata, handles)
global Generator
Generator.ChargeEff=str2double(get(hObject,'string'))/100;

% --- Executes during object creation, after setting all properties.
function editChargeEff_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editDischargeEff_Callback(hObject, eventdata, handles)
global Generator
Generator.DischargeEff=str2double(get(hObject,'string'))/100;

% --- Executes during object creation, after setting all properties.
function editDischargeEff_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editLoss_Callback(hObject, eventdata, handles)
% hObject    handle to editLoss (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Generator
Generator.VariableStruct.SelfDischarge = str2double(get(hObject,'string'))/100/24;

% --- Executes during object creation, after setting all properties.
function editLoss_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editLoss (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function fillRatesStatic(handles)
global Plant Generator

demand = Plant.Data.Demand;
ChillSizeTons = 0;
kWh2ton = 1/3.51685;
for i = 1:1:length(Plant.Generator)
    if strcmp(Plant.Generator(i).Type,'Chiller')
        ChillSizeTons  = ChillSizeTons +Plant.Generator(i).Size;
    end
end
ChillSizekWh = ChillSizeTons/kWh2ton;
% calculate the TES fill and discharge rates using demand & chiller capacity
FillRatePerc = 100/(Generator.Size/ChillSizekWh); % 100% divided by # of hours to fill
FillRate = Generator.VariableStruct.SizeLiter*FillRatePerc/60; % fill rate in liter/min
if strcmp(Generator.VariableStruct.EnStoreType, 'ColdTES')
    MaxDemandC = max(demand.C);
    DischRatePerc = (1000*MaxDemandC)/Generator.Size*100;
    DischRate = Generator.VariableStruct.SizeLiter*DischRatePerc/60; % disch rate in liter/min
elseif strcmp(Generator.VariableStruct.EnStoreType, 'HotTES')
    MaxDemandH = max(demand.H); 
    DischRatePerc = (1000*MaxDemandH)/Generator.Size*100;
    DischRate = Generator.VariableStruct.SizeLiter*DischRatePerc/60; % disch rate in liter/min
end
set(handles.textFillRate,'string',FillRate)
set(handles.textFillRatePerc,'string',FillRatePerc)
set(handles.textDischRate,'string',DischRate)
set(handles.textDischRatePerc,'string',DischRatePerc)


% --- Executes on button press in pushbuttonComPorts.
function pushbuttonComPorts_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonComPorts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global EditCommHandle
CommunicationPorts();
waitfor(EditCommHandle)
