function varargout = SetupHvac(varargin)
% SETUPHVAC MATLAB code for SetupHvac.fig
%      SETUPHVAC, by itself, creaHVAC a new SETUPHVAC or raises the existing
%      singleton*.
%
%      H = SETUPHVAC returns the handle to a new SETUPHVAC or the handle to
%      the existing singleton*.
%
%      SETUPHVAC('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SETUPHVAC.M with the given input arguments.
%
%      SETUPHVAC('Property','Value',...) creaHVAC a new SETUPHVAC or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SetupHvac_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SetupHvac_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SetupHvac

% Last Modified by GUIDE v2.5 21-Mar-2016 13:47:41

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SetupHvac_OpeningFcn, ...
                   'gui_OutputFcn',  @SetupHvac_OutputFcn, ...
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


% --- ExecuHVAC just before SetupHvac is made visible.
function SetupHvac_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SetupHvac (see VARARGIN)

% Choose default command line output for SetupHvac
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

global Plant SYSINDEX figHandle Generator
figHandle = handles.figure1;
Generator = Plant.Generator(SYSINDEX);
HVAC=Plant.Generator(SYSINDEX).VariableStruct;

set(handles.editName,'string',Generator.Name)
% set(handles.textType,'string',COMPONENT.Type)
set(handles.editSize,'string',Generator.Size)
if strcmp(HVAC.EnStoreType, 'AC')
    set(handles.uipanelMode,'SelectedObject', handles.radiobuttonAC)
elseif strcmp(HVAC.EnStoreType, 'HV')
    set(handles.uipanelMode,'SelectedObject', handles.radiobuttonHV)
elseif strcmp(HVAC.EnStoreType, 'HVAC')
    set(handles.uipanelMode, 'SelectedObject', handles.radiobuttonHVAC)
end

%HVAC specific
%set(handles.editSize,'string',HVAC.SizeLiter)
set(handles.editTcold,'string',HVAC.Tcold)%this is the coldest you will allow your building to get
set(handles.editThot,'string',HVAC.Thot)%this is the hottest you will allow your building to get
set(handles.editLoss,'string',HVAC.SelfDischarge*100*24)%thermal energy lost per day
set(handles.footageSize, 'string', HVAC.footageSize)
%set(handles.editFillRate,'string',HVAC.FillRate)%this is equal to the output of the chiller or heater
%set(handles.editFillRatePerc,'string',HVAC.FillRatePerc)
%set(handles.editDischRate,'string',HVAC.DischRate)%this rate limit is set to infinity and cannot be changed because there the 'storage' is immediately sent to the demand
%set(handles.editDischRatePerc,'string',HVAC.DischRatePerc)
%set(handles.editDischargeEff,'string',HVAC.DischargeEff*100)
%set(handles.editChargeEff,'string',HVAC.ChargeEff*100)

fillRatesStatic(handles)


% --- Outputs from this function are returned to the command line.
function varargout = SetupHvac_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


%the size that the user inputs should be the size of the area that is
%temperature controlled with the smart HVAC system.
function footageSize_Callback(hObject, eventdata, handles)
% hObject    handle to footageSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Generator
Litersperft3 = 28.3618; %the user input should be in square footage. multiply that by 10ft tall for cubic feet, then use ft3toLiter to find volume in liters
kWh2Liter = 3600/1.005; %3600 seconds / 1.005 kJ/kg*K
densityair = 1.194; %kg/m^3 density of air at ~23C
Cpair = 1.005; %constant pressure specific heat of air at ~20C in kJ/kgK
Generator.VariableStruct.footageSize=str2double(get(hObject,'string'));
%size is mCPdeltaT
%this is the size in kWh. footagesize*10=cubicfeet; Litersperft3*.001 turns
%cubic feet to m^3;densityair takes it to mass; Cpair takes it to kJ/kgK;
%deltaT takes it to kJ; dividing by 3600 converts to kWh from kJ
Generator.Size=Generator.VariableStruct.footageSize*10*(Litersperft3*.001)*densityair*Cpair*(Generator.VariableStruct.Thot-Generator.VariableStruct.Tcold)/3600;%this is the size in kWh
Generator.VariableStruct.SizeLiter = Generator.Size*kWh2Liter/(Generator.VariableStruct.Thot-Generator.VariableStruct.Tcold);%this is the size in liters
set(handles.editSize,'string',Generator.Size)


% --- ExecuHVAC during object creation, after setting all properties.
function footageSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to footageSize (see GCBO)
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


% --- ExecuHVAC during object creation, after setting all properties.
function editName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- ExecuHVAC on button press in pushbuttonSaveOnlyToPlant.
function pushbuttonSaveOnlyToPlant_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonSaveOnlyToPlant (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Plant SYSINDEX Generator
Plant.Generator(SYSINDEX) = Generator;
clear global Generator
close(gcf)

% --- ExecuHVAC on button press in pushbuttonCancel.
function pushbuttonCancel_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonCancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close(gcf)

% --- ExecuHVAC on button press in pushbuttonSaveAs.
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
savedir=fullfile(Model_dir,'component library','HVAC',strcat(Generator.Name,'.mat'));
[f,p]=uiputfile(savedir,'Save As Energy Storage Component');
if f==0;return;end
save([p f],'component')
pushbuttonSaveOnlyToPlant_Callback(hObject, eventdata, handles)


%you can also input the size in liters and get out the size in kWh
function editSize_Callback(hObject, eventdata, handles)
% hObject    handle to editSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Generator
Generator.Size=str2double(get(hObject,'string'));
kWh2Liter = 3600/1.005; %3600 seconds / 1.005 kJ/kg*K
Litersperft3 = 28.3618;
Cpair = 1.005;%kJ/kgK
densityair = 1.194; %kg/m3
Generator.VariableStruct.SizeLiter = Generator.Size*kWh2Liter/(Generator.VariableStruct.Thot-Generator.VariableStruct.Tcold);
%assuming 10 feet tall, convert storage size to size of building
%size*3600 turns kWh to kJ. divide by cpdeltaT for mass
%divide by density for volume. multiply by 1000 to go from m3 to liters. convert to ft and divide by 10 feet tall.
Generator.VariableStruct.footageSize = Generator.Size*3600/(Cpair*(Generator.VariableStruct.Thot-Generator.VariableStruct.Tcold))/densityair*1000/Litersperft3;
set(handles.footageSize,'string',Generator.VariableStruct.footageSize)


% --- ExecuHVAC during object creation, after setting all properties.
function editSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editSize (see GCBO)
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
kWh2Liter = 3600/1.005; %3600 seconds / 1.005 kJ/kg*K 
Generator.Size = Generator.VariableStruct.SizeLoadShift/kWh2Liter*(Generator.VariableStruct.Thot-Generator.VariableStruct.Tcold);
set(handles.editSize,'string',Generator.Size)


% --- ExecuHVAC during object creation, after setting all properties.
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
kWh2Liter = 3600/1.005; %3600 seconds / 1.005 kJ/kg*K 
Generator.Size = Generator.VariableStruct.SizeLiter/kWh2Liter*(Generator.VariableStruct.Thot-Generator.VariableStruct.Tcold);
set(handles.editSize,'string',Generator.Size)


% --- ExecuHVAC during object creation, after setting all properties.
function editThot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editThot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function editLoss_Callback(hObject, eventdata, handles)
% hObject    handle to editLoss (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Generator
Generator.VariableStruct.SelfDischarge = str2double(get(hObject,'string'))/100/24;

% --- ExecuHVAC during object creation, after setting all properties.
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
HeatSizeTons = 0;
kWh2ton = 1/3.51685;
for i = 1:1:length(Plant.Generator)
    if strcmp(Plant.Generator(i).Type,'Chiller')
        ChillSizeTons  = ChillSizeTons +Plant.Generator(i).Size;
    end
    if strcmp(Plant.Generator(i).Type, 'Heater')
        HeatSizeTons = HeatSizeTons + Plant.Generator(i).Size;
    end
end
ChillSizekWh = ChillSizeTons/kWh2ton;
HeatSizekWh = HeatSizeTons/kWh2ton;
% calculate the HVAC fill rate of HVAC using demand & chiller capacity
% the fill rate percent is solely dependent on how quickly the heater or
% chiller can fill it.
%FillRatePerc = max(100*(ChillSizekWh/Generator.Size), 100*(HeatSizekWh/Generator.Size)); % 100% divided by # of hours to fill
FillRatePerc =100;
FillRate = inf;
%FillRate = Generator.VariableStruct.SizeLiter*FillRatePerc/60; % fill rate in liter/min. for HVAC this is the entire size times the percent fill rate
Generator.VariableStruct.FillRatePerc = FillRatePerc;
Generator.VariableStruct.FillRate = FillRate;
if strcmp(Generator.VariableStruct.EnStoreType, 'AC')
    MaxDemandC = max(demand.C);
    DischRatePerc = (1000*MaxDemandC)/Generator.Size*100;
    DischRate = inf; % disch rate in liter/min
elseif strcmp(Generator.VariableStruct.EnStoreType, 'HV')
    MaxDemandH = max(demand.H); 
    DischRatePerc = (1000*MaxDemandH)/Generator.Size*100;
    DischRate = inf; % disch rate in liter/min
elseif strcmp(Generator.VariableStruct.EnStoreType, 'HVAC')
    MaxDemandH = max(demand.H);
    MaxDemandC = max(demand.C);
    DischRatePerc = (1000*max(MaxDemandH, MaxDemandC))/Generator.Size*100;
    DischRate = inf;
end
set(handles.textFillRate,'string',FillRate)
set(handles.textFillRatePerc,'string',FillRatePerc)
% set(handles.textDischRate,'string',DischRate)
% set(handles.textDischRatePerc,'string',DischRatePerc)


% --- ExecuHVAC on button press in pushbuttonComPorts.
function pushbuttonComPorts_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonComPorts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global EditCommHandle
CommunicationPorts();
waitfor(EditCommHandle)


% --- Executes when selected object is changed in uipanelMode.
function uipanelMode_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanelMode 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

