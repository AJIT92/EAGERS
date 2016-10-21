function varargout = EAGERS(varargin)
% EAGERS MATLAB code for EAGERS.fig
%      EAGERS, by itself, creates a new EAGERS or raises the existing
%      singleton*.
%
%      H = EAGERS returns the handle to a new EAGERS or the handle to
%      the existing singleton*.
%
%      EAGERS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EAGERS.M with the given input arguments.
%
%      EAGERS('Property','Value',...) creates a new EAGERS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before EAGERS_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to EAGERS_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help EAGERS

% Last Modified by GUIDE v2.5 13-Sep-2016 21:39:03
    
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @EAGERS_OpeningFcn, ...
                   'gui_OutputFcn',  @EAGERS_OutputFcn, ...
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


function EAGERS_OpeningFcn(hObject, eventdata, handles, varargin)
% Choose default command line output for EAGERS
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

global Model_dir Plant BUILDING figHandle HistProf Holidays
figHandle = [];
Plant = [];
files=dir(fullfile(Model_dir, 'Plant','*.mat'));
list=strrep({files.name},'.mat','');
if isempty(Plant)
    K = menu('Opening Dialog','Open existing Project','Start New Project');
    if K ==1
        answer = char(list(listdlg('PromptString','Select project','SelectionMode','single','ListString',list)));
        load(fullfile(Model_dir,'Plant',answer))
    elseif K==2
       %create project
       Name = char(inputdlg('Specify the project name','Project Name', 1,{'MicroGrid_01'}));
       BUILDING = [];
       hCreate = dialog('Visible','off');
       CreateNewProject(Name,hCreate)
       waitfor(hCreate)
    end 
end    
set(gcf,'Name','EAGERS 2016.9.13')

movegui(gcf,'center');

DataNames1 = cellstr(strcat('Demand.',fields(Plant.Data.Demand)));
set(handles.popupmenuData,'string',[cellstr('Temperature');DataNames1],'value',1)

names = {};
for i= 1:1:length(Plant.Generator)
    names(end+1) = {Plant.Generator(i).Name};
end
set(handles.popupmenuGenerators,'string',names,'value',1)

value=strmatch(Plant.Name,list,'exact');
set(handles.popupmenuPlant,'string',list,'value',value)
popupmenuPlant_Callback(handles.popupmenuPlant,eventdata,handles,1)

set(handles.textBaseline,'string','Not Yet Calculated')
set(handles.textDispatch,'string','Not Yet Calculated')

days = round(Plant.Data.Timestamp(end) - Plant.Data.Timestamp(1));
set(handles.sliderZoom1,'Min',1,'Max',4,'Value',1,'SliderStep',[1/3,1/3])
set(handles.sliderDate1,'Min',1,'Max',2,'Value',1,'SliderStep',[1/(days-1),1/(days-1)])
sliderZoom1_Callback(handles.sliderZoom1, eventdata, handles)
if isfield(Plant,'Dispatch')
    if isfield(Plant.Dispatch,'Historical')
        set(handles.sliderDate1,'Min',1,'Max',2,'Value',1,'SliderStep',[1/(days-1),1/(days-1)])
    else
        days = round(Plant.Dispatch.Dispatch.Timestamp(end) - Plant.Dispatch.Dispatch.Timestamp(1));
        set(handles.sliderDate2,'Min',1,'Max',2,'Value',1,'SliderStep',[1/days,1/days])
    end
    popupmenuGenerators_Callback(handles.popupmenuGenerators,eventdata,handles)
end


if isfield(Plant.Data,'Holidays')
    Holidays = Plant.Data.Holidays;
else Holidays =[];
end
if isfield(Plant.Data,'HistProf')
    HistProf = Plant.Data.HistProf;
else HistProf =[];
end
if ~isfield(Plant.optimoptions,'Buffer')
    Plant.optimoptionsBuffer = 20; % percentage for buffer on storage
end

function varargout = EAGERS_OutputFcn(hObject, eventdata, handles)
varargout{1} = handles.output;


function Plot2Axes(handles)
global Plant
val=get(handles.popupmenuData,'value');
DataNames=get(handles.popupmenuData,'string');
currentstr=DataNames{val};
[PlotIndex, forecast] = PlotWindow(handles, Plant.Data.Timestamp, currentstr);
TimeVec = Plant.Data.Timestamp(PlotIndex(1):PlotIndex(2))';
if strcmp(currentstr,'Temperature')
    currentstr ='T';
    Y = Plant.Data.Temperature(PlotIndex(1):PlotIndex(2))';
elseif strncmp(currentstr,'Demand',6)
    Y = Plant.Data.Demand.(currentstr(end))(PlotIndex(1):PlotIndex(2));
    currentstr=currentstr(end);
end
axes(handles.axes1)
PlotData(TimeVec,Y,forecast,currentstr,Plant.Data,[],[],[])

function [PlotIndex, forecast] = PlotWindow(handles, timestamp,currentstr)
a = datevec(timestamp(1));
b = datevec(timestamp(end));
D1 = datenum([a(1) a(2) a(3)]);
months = max(1,12*(b(1)-a(1))+b(2)-a(2));
years = max(1,b(1)-a(1));
days = round(timestamp(end)-timestamp(1));
weeks = ceil(days/7);
Zoom = max(1,round(get(handles.sliderZoom1,'Value')));
if Zoom == 1
    day = 1+round((get(handles.sliderDate1,'Value')-1)*(days-1));
    endday = day;
elseif Zoom ==2
    day = 1+7*round((get(handles.sliderDate1,'Value')-1)*(weeks-1));
    endday = day+6;
elseif Zoom ==3
    month = 1+ round((get(handles.sliderDate1,'Value')-1)*(months-1));
    day = datenum(a(1),month,1)-datenum(a(1),1,1)+1;
    endday = day+(datenum(a(1),month+1,1)-datenum(a(1),month,1)-1);
elseif Zoom ==4
    year = a(1)+ round((get(handles.sliderDate1,'Value')-1)*(years-1));
    day = 1+round(datenum(year,1,1)-round(timestamp(1)));
    endday = round(datenum(year+1,1,1)-round(timestamp(1)));
end
PlotIndex(1) = max(1,nnz(timestamp<=(D1+(day-1))));
PlotIndex(2) = nnz(timestamp<=(D1+endday));
if PlotIndex(2)>length(timestamp)
    set(handles.sliderZoom,'Value',Zoom-1)
    PlotWindow(handles, timestamp, currentstr)
end
if (strcmp(currentstr,'Temperature') || strncmp(currentstr,'Demand',6)) && (endday-day)<7 && get(handles.checkboxForecast,'Value')
    forecast = 1;
else forecast =0;
end

function popupmenuGenerators_Callback(hObject, eventdata, handles)
global Plant
val=get(handles.popupmenuGenerators,'value');
Timestamp = Plant.Data.Timestamp;
days = round(Timestamp(end)-Timestamp(1));
day = 1+round((get(handles.sliderDate2,'Value')-1)*(days-1));
Timeplot = Timestamp(1)+(day-1);%date number of plot start
axes(handles.axes2)
HistDisp=[];
RunDisp =[];
LEG ={};
X=[];
Y =[];
currentstr = [];
if isfield(Plant,'Dispatch')
    if isfield(Plant.Dispatch,'Dispatch') 
        Timestamp = Plant.Dispatch.Dispatch.Timestamp;
        Ts2 = (Timestamp(3)-Timestamp(2))*24;%time of 1 step in hours
        if (Plant.Dispatch.Dispatch.Timestamp(1)<=Timeplot && Plant.Dispatch.Dispatch.Timestamp(end)>Timeplot)
            Xs = nnz(Timestamp<=Timeplot);
            XFs = nnz(Timestamp<=Timeplot+1);
            X = Timestamp(Xs:XFs);
            Y = 0*X;
        elseif ~isfield(Plant.Dispatch,'Historical')
            Xs = 1;
            XFs = min(length(Timestamp),nnz(Timestamp<=(Timestamp(1)+1)));
            X = Timestamp(Xs:XFs);
            Y = 0*X;
        end
        if ~isempty(Y)
            currentstr = 'Disp';
            LEG={};
            if (~isempty(strfind(Plant.Generator(val).Type,'Storage'))||~isempty(strfind(Plant.Generator(val).Type,'Hvac'))) %might neet to add gentype 12 for thermal energy storage from Hvac
                Y(2:end,1)=zeros(XFs-Xs,1);%there is no storage in the baseline %(Plant.Dispatch.Baseline.GeneratorState(Xs:XFs-1,val)-Plant.Dispatch.Baseline.GeneratorState(Xs+1:XFs,val))/Ts2;%energy storage power output in kW
%               Y(2:end,2)=(Plant.Dispatch.Predicted(2).GenDisp(Xs:XFs-1,val)-Plant.Dispatch.Predicted(2).GenDisp(Xs+1:XFs,val))/Ts2;%energy storage power output in kW
                if get(handles.checkboxEnergyStorageState,'Value')==1
                    currentstr = 'SOC';
                end
            else
            %% Baseline is not currently calculated
%             if length(Plant.Dispatch.Baseline.Timestamp)<length(X) %this occurs during variable timesteps or short horizons
%                 Y(:,1)=interp1([Plant.Dispatch.Baseline.Timestamp,X(end)],[Plant.Dispatch.Baseline.GeneratorState(:,val);Plant.Dispatch.Baseline.GeneratorState(1,val)],X);
%                 %Y(:,2)=interp1([Plant.Dispatch.Predicted.Timestamp,X(end)],[Plant.Dispatch.Predicted.GeneratorState(:,val);Plant.Dispatch.Predicted.GeneratorState(1,val)],X);
%             else Y(:,1)=Plant.Dispatch.Baseline.GeneratorState(Xs:XFs,val);
%             end
%             Y(:,2)=Plant.Dispatch.Predicted(2).GenDisp(Xs:XFs,val);
%             LEG={'Baseline','Predicted'};
            end
        end
    end
end

if isfield(Plant,'Dispatch')
    if isfield(Plant.Dispatch,'Historical')
        Timestamp = Plant.Data.Timestamp;
        Ts = (Timestamp(2)-Timestamp(1))*24;%time of 1 step in hours
        Xs = nnz(Timestamp<=Timeplot);
        XFs = nnz(Timestamp<=Timeplot+1);
        HistDisp.X = Timestamp(Xs:XFs);
        if (~isempty(strfind(Plant.Generator(val).Type,'Storage'))||~isempty(strfind(Plant.Generator(val).Type,'Hvac'))) %energy storage
            if get(handles.checkboxEnergyStorageState,'Value')
                HistDisp.Y = [0*HistDisp.X 0*HistDisp.X];
                HistDisp.Y(:,2) = Plant.Dispatch.Historical(Xs:XFs,val);
            else HistDisp.Y = 0*HistDisp.X;
            end
            for t = max(Xs,2):1:XFs %avoid calculating @ t= 1
                HistDisp.Y(t,1) = (Plant.Dispatch.Historical(t,val)-Plant.Dispatch.Historical(t-1,val))/Ts; %energy storage power output in kW
            end
            LEG(end+1:end+2) = {'Historical (kW)','Historical SOC (kWh)'};
        else HistDisp.Y = Plant.Dispatch.Historical(Xs:XFs,val);%not energy storage
            LEG(end+1)={'Historical'};
        end
    end
end

if isfield(Plant,'Dispatch')
    if isfield(Plant.Dispatch,'RunData')
        Timestamp = Plant.Dispatch.RunData.Timestamp;
        Ts = (Timestamp(2)-Timestamp(1))*24;%time of 1 step in hours
        if Timestamp(1)<=Timeplot && Timestamp(end)>Timeplot
            Xs = nnz(Timestamp<=Timeplot);
            XFs = nnz(Timestamp<=Timeplot+1);
        else
            Xs = 1;
            XFs = min(length(Timestamp),nnz(Timestamp<=(Timestamp(1)+1)));
        end
        RunDisp.X = Timestamp(Xs:XFs);
        if (~isempty(strfind(Plant.Generator(val).Type,'Storage'))||~isempty(strfind(Plant.Generator(val).Type,'Hvac'))) %energy storage
            if get(handles.checkboxEnergyStorageState,'Value')
                RunDisp.Y = [0*RunDisp.X 0*RunDisp.X];
                for t = max(Xs,2):1:XFs %avoid calculating @ t= 1
                    RunDisp.Y(t,1) = (Plant.Dispatch.RunData.GeneratorState(t-1,val)-Plant.Dispatch.RunData.GeneratorState(t,val))/Ts; %energy storage power output in kW
                end
                RunDisp.Y(:,2) = Plant.Dispatch.RunData.GeneratorState(Xs:XFs,val);
                LEG = {'Dispatch (kW)'; 'Dispatch SOC (kWh)'};
            else RunDisp.Y = 0*RunDisp.X;
                for t = max(Xs,2):1:XFs %avoid calculating @ t= 1
                    RunDisp.Y(t,1) = (Plant.Dispatch.RunData.GeneratorState(t-1,val)-Plant.Dispatch.RunData.GeneratorState(t,val))/Ts; %energy storage power output in kW
                end
                LEG(end+1) = {'Dispatch (kW)'};
            end
        else RunDisp.Y = Plant.Dispatch.RunData.GeneratorState(Xs:XFs,val);%not energy storage
            LEG(end+1)={'Dispatch'};
        end
    end
end
if ~isempty(X)||~isempty(RunDisp)||~isempty(HistDisp)
    PlotData(X,Y,0,currentstr,Plant.Data,HistDisp,RunDisp,[])
    legend(LEG)   
end

function sliderZoom1_Callback(hObject, eventdata, handles)
global Plant
Zoom = round(get(hObject,'Value'));
if Zoom <=2
     set(handles.checkboxForecast,'enable','on')
else set(handles.checkboxForecast,'enable','off')
end
a = datevec(Plant.Data.Timestamp(1));
b = datevec(Plant.Data.Timestamp(end));
months = max(1,12*(b(1)-a(1))+b(2)-a(2));
years = max(1,b(1)-a(1));
days = Plant.Data.Timestamp(end)-Plant.Data.Timestamp(1)+1;

if Zoom ==1
    set(handles.sliderDate1,'SliderStep',[1/days,10/days])
elseif Zoom == 2
    set(handles.sliderDate1,'SliderStep',[7/days,70/days])
elseif Zoom == 3
    set(handles.sliderDate1,'SliderStep',[1/months,1/months])
else set(handles.sliderDate1,'SliderStep',[1/years,1/years])
end
Plot2Axes(handles)

function popupmenuPlant_Callback(hObject, eventdata, handles, skipreload)
global Plant Model_dir Holidays HistProf
PlantList=get(hObject,'string');
PlantVal=get(hObject,'value');
if nargin ==3 || skipreload == 0
    projdir=fullfile(Model_dir, 'Plant');
    projName=fullfile(projdir,PlantList{PlantVal});
    load (projName)
end
if isfield(Plant.Data,'Holidays')
    Holidays = Plant.Data.Holidays;
else Holidays =[];
end
if isfield(Plant.Data,'HistProf')
    HistProf = Plant.Data.HistProf;
else HistProf =[];
end
listbox_MakeList(hObject, eventdata, handles)
days = round(Plant.Data.Timestamp(end) - Plant.Data.Timestamp(1));
set(handles.sliderZoom1,'Min',1,'Max',4,'Value',1,'SliderStep',[1/3,1/3])
set(handles.sliderDate1,'Min',1,'Max',2,'Value',1,'SliderStep',[1/(days-1),1/(days-1)])
sliderZoom1_Callback(handles.sliderZoom1, eventdata, handles)
if isfield(Plant,'Dispatch')
    if isfield(Plant.Dispatch,'Historical')
        Timestamp = Plant.Data.Timestamp;
    else Timestamp = Plant.Dispatch.Dispatch.Timestamp;
    end
    days = round(Timestamp(end) - Timestamp(1));
    set(handles.sliderDate2,'Min',1,'Max',2,'Value',1,'SliderStep',[1/days,1/days])
    names = {};
    for i = 1:1:length(Plant.Generator)
        names(end+1) = {Plant.Generator(i).Name};
    end
    set(handles.popupmenuGenerators,'string',names,'value',1)
    popupmenuGenerators_Callback(handles.popupmenuGenerators,eventdata,handles)
else%even if there is no dispatch set dropdown correctly
    names = {};
    for i= 1:1:length(Plant.Generator)
        names(end+1) = {Plant.Generator(i).Name};
    end
    set(handles.popupmenuGenerators,'string',names,'value',1)
end
if ~isfield(Plant.optimoptions,'Buffer')
    Plant.optimoptionsBuffer = 20; % percentage for buffer on storage
end
DataNames1 = cellstr(strcat('Demand.',fields(Plant.Data.Demand)));
set(handles.popupmenuData,'string',[cellstr('Temperature');DataNames1],'value',1)
set(handles.textBaseline,'string','Not Yet Calculated')
set(handles.textDispatch,'string','Not Yet Calculated')

function pushbuttonPlant_Callback(hObject, eventdata, handles)
global Plant Model_dir
[f,p]=uiputfile(fullfile(Model_dir,'Plant','PlantNew.mat'),'Save Plant As...');
if f==0; return; end
Plant.id=str2double(datestr(now,'yymmddHHMMSS'));
Plant.Name=strrep(f,'.mat','');
save([p,f],'Plant')
listbox_MakeList(hObject, eventdata, handles)
files=dir(fullfile(Model_dir, 'Plant','*.mat'));
list=strrep({files.name},'.mat','');
value=nonzeros((1:length(list)).*strcmp(Plant.Name,list));
set(handles.popupmenuPlant,'string',list,'value',value)

function listbox_MakeList(hObject, eventdata, handles)
global Plant
Gen = [];
chpGen = Gen;
Chill = Gen;
Heat = Gen;
Renew = Gen;
Storage = Gen;
Util = Gen;
Hvac = Gen;
str={};
for i=1:length(Plant.Generator)
    if strcmp(Plant.Generator(i).Type,'Utility')
        Util(end+1) = i;
    elseif strcmp(Plant.Generator(i).Type,'CHP Generator')
        chpGen(end+1) = i;
    elseif strcmp(Plant.Generator(i).Type,'Electric Generator')
        Gen(end+1) = i;
    elseif strcmp(Plant.Generator(i).Type,'Chiller')
        Chill(end+1) = i;
    elseif strcmp(Plant.Generator(i).Type,'Heater') || strcmp(Plant.Generator(i).Type,'Boiler')
        Heat(end+1) = i;
    elseif strcmp(Plant.Generator(i).Type,'Solar')
        Renew(end+1) = i;
    elseif strcmp(Plant.Generator(i).Type,'Wind')
        Renew(end+1) = i;
    elseif strcmp(Plant.Generator(i).Type,'Thermal Storage') || strcmp(Plant.Generator(i).Type,'Electric Storage')
        Storage(end+1) = i;
    elseif strcmp(Plant.Generator(i).Type,'Hvac')
        Hvac(end+1) = i;
    end
end
if nnz(Util)>0
    str(end+1)={'Utilities:'};
end
for i=1:1:nnz(Util)
    str(end+1) = {['      ' Plant.Generator(Util(i)).Name]};
    if Plant.Generator(Util(i)).Enabled ==0
      str{end}=strcat(str{end},'--Offline');
    end
end
if nnz(chpGen)>0
    str(end+1)={'CHP Generators:'};
end
for i=1:1:nnz(chpGen)
    str(end+1) = {['      ' Plant.Generator(chpGen(i)).Name]};
    if Plant.Generator(chpGen(i)).Enabled ==0
      str{end}=strcat(str{end},'--Offline');
    end
end
if nnz(Gen)>0
    str(end+1)={'Electric Generators:'};
end
for i=1:1:nnz(Gen)
    str(end+1) = {['      ' Plant.Generator(Gen(i)).Name]};
    if Plant.Generator(Gen(i)).Enabled ==0
      str{end}=strcat(str{end},'--Offline');
    end
end
if nnz(Heat)>0
    str(end+1)={'Heaters/Boilers:'};
end
for i=1:1:nnz(Heat)
    str(end+1) = {['      ' Plant.Generator(Heat(i)).Name]};
    if Plant.Generator(Heat(i)).Enabled ==0
      str{end}=strcat(str{end},'--Offline');
    end
end
if nnz(Chill)>0
    str(end+1)={'Chillers:'};
end
for i=1:1:nnz(Chill)
    str(end+1) = {['      ' Plant.Generator(Chill(i)).Name]};
    if Plant.Generator(Chill(i)).Enabled ==0
      str{end}=strcat(str{end},'--Offline');
    end
end
if nnz(Renew)>0
    str(end+1)={'Renewables:'};
end
for i=1:1:nnz(Renew)
    str(end+1) = {['      ' Plant.Generator(Renew(i)).Name]};
    if Plant.Generator(Renew(i)).Enabled ==0
      str{end}=strcat(str{end},'--Offline');
    end
end
if nnz(Storage)>0
    str(end+1)={'Energy Storage:'};
end
for i=1:length(Storage)
    str(end+1) = {['      ' Plant.Generator(Storage(i)).Name]};
    if Plant.Generator(Storage(i)).Enabled ==0
      str{end}=strcat(str{end},'--Offline');
    end
end
if nnz(Hvac)>0
    str(end+1)={'Smart HVAC:'};
end
for i=1:length(Hvac)
    str(end+1) = {['      ' Plant.Generator(Hvac(i)).Name]};
    if Plant.Generator(Hvac(i)).Enabled ==0
      str{end}=strcat(str{end},'--Offline');
    end
end
set(handles.listboxSetup,'string',str,'value',1)

function listboxSetup_Callback(hObject, eventdata, handles)
list=get(handles.listboxSetup,'string');
val=get(handles.listboxSetup,'value');
compSel=list{val};
if strncmp(compSel,'    ',4)
   set(handles.pushbuttonEdit,'enable','on')
   set(handles.pushbuttonRemove,'enable','on')
   set(handles.pushbuttonAdd,'enable','on')
   set(handles.pushbuttonEnable,'enable','on')
   set(handles.pushbuttonDisable,'enable','on')
else
   set(handles.pushbuttonEdit,'enable','off')
   set(handles.pushbuttonRemove,'enable','off')
   set(handles.pushbuttonAdd,'enable','off')
   set(handles.pushbuttonEnable,'enable','off')
   set(handles.pushbuttonDisable,'enable','off')
end  

function listboxSetup_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over listboxSetup.
function listboxSetup_ButtonDownFcn(hObject, eventdata, handles)
global Plant SYSINDEX figHandle
list=get(handles.listboxSetup,'string');
listnum=get(handles.listboxSetup,'value');
selection=list{listnum};
if strncmp(selection,'      ',6)
    selection = strrep(selection,'      ','');
    GenNames = {};
    for i = 1:1:length(Plant.Generator)
        GenNames{i} = char(Plant.Generator(i).Name);
    end
    SYSINDEX = find(strcmp(selection,GenNames),1);
end
if max(strcmp(Plant.Generator(SYSINDEX).Type,{'Electric Generator';'CHP Generator';'Chiller';'Heater';'Boiler';}))==1
    SetupGenerator
elseif strcmp(Plant.Generator(SYSINDEX).Type,'Thermal Storage')
    SetupThermalStorage
elseif strcmp(Plant.Generator(SYSINDEX).Type,'Electric Storage')
    SetupBatteryStorage
elseif strcmp(Plant.Generator(SYSINDEX).Type,'Utility')
    if strcmp(Plant.Generator(SYSINDEX).Source,'Electricity')
        SetupUtility('Electric')
    else
        SetupUtility('Gas')
    end
elseif strcmp(Plant.Generator(SYSINDEX).Type,'Solar')
    SetupSolar
elseif strcmp(Plant.Generator(SYSINDEX).Type,'Wind')
    %    SetupWind
elseif strcmp(Plant.Generator(SYSINDEX).Type,'Hvac')
    SetupHvac
end
waitfor(figHandle)
% if isfield(Plant.Generator,'OpMatA')%if you have eddited the generator after a run then 
%     %OpMatA will already exist so clear it, otherwise, it won't exist so
%     %don't create it
%     Plant.Generator(SYSINDEX).OpMatA = {};
% end
listbox_MakeList(hObject, eventdata, handles)


function pushbuttonAdd_Callback(hObject, eventdata, handles)
global Model_dir Plant
list={'Utility'
    'Electric Generator'
    'CHP Generator'
    'Chiller'
    'Heater'
    'Boiler'
    'Wind'
    'Solar'
    'Thermal Storage'
    'Electric Storage'
    'HVAC'};

[s,v] = listdlg('PromptString','Select System Category', 'SelectionMode','single','ListString',list);
compdir=fullfile(Model_dir, 'component library', char(list(s))); 
files=dir(fullfile(compdir,'*.mat'));
list2=strrep({files.name},'.mat','');
[s2,v] = listdlg('PromptString','Select System', 'SelectionMode','single','ListString',list2);
filename = fullfile(Model_dir, 'component library',char(list(s)),strcat(char(list2(s2)),'.mat'));
load(filename)
if isfield(component,'SumStartMonth') %this is if you are adding a new electric utility with a rate table
    comp.Type = component.Type;
    comp.Name = component.Name;
    comp.Source = 'Electricity';
    comp.Output = [];
    comp.Output.E = 1;
    comp.Size = 1;
    comp.Enabled = 1;
    comp.VariableStruct = [];
    comp.VariableStruct = component;
    component = comp;
end
if isfield(Plant.Generator(1),'OpMatA')
    component.OpMatA = {};%if you are adding a new component add these fields
end
if isfield(Plant.Generator(1),'OpMatB')
    component.OpMatB = {};
end
Plant.Generator(end+1) = component;
popupmenuPlant_Callback(handles.popupmenuPlant,eventdata,handles,1);

function pushbuttonEdit_Callback(hObject, eventdata, handles)
listboxSetup_ButtonDownFcn(hObject, eventdata, handles)

function pushbuttonRemove_Callback(hObject, eventdata, handles)
global Plant
names = {};
for i = 1:1:length(Plant.Generator)
    names(end+1) = cellstr(Plant.Generator(i).Name);
end
list=get(handles.listboxSetup,'string');
val=get(handles.listboxSetup,'value');
componentSelected=list(val);
componentSelected = strrep(char(componentSelected),'      ','');
answer=questdlg({'Are you sure you want to remove:',componentSelected});
if strcmp(answer,'Yes')
    GenNum = find(strcmp(componentSelected,names));
    Gen2 = Plant.Generator(1:GenNum-1);
    Gen2(end+1:length(Plant.Generator)-1) = Plant.Generator(GenNum+1:end);
    Plant.Generator = Gen2;
    popupmenuPlant_Callback(handles.popupmenuPlant,eventdata,handles,1);
end

% --- Executes on button press in pushbuttonDisable.
function pushbuttonDisable_Callback(hObject, eventdata, handles)
global Plant
names = {};
for i = 1:1:length(Plant.Generator)
    names(end+1) = cellstr(Plant.Generator(i).Name);
end
list=get(handles.listboxSetup,'string');
listnum=get(handles.listboxSetup,'value');
selection=list{listnum};
selection = strrep(selection,'      ','');
GenNum = find(strcmp(selection,names));
Plant.Generator(GenNum).Enabled = 0;

listbox_MakeList(hObject, eventdata, handles)

% --- Executes on button press in pushbuttonEnable.
function pushbuttonEnable_Callback(hObject, eventdata, handles)
global Plant
names = {};
for i = 1:1:length(Plant.Generator)
    names(end+1) = cellstr(Plant.Generator(i).Name);
end
list=get(handles.listboxSetup,'string');
listnum=get(handles.listboxSetup,'value');
selection=list{listnum};
selection = strrep(selection,'      ','');
selection = strrep(selection,'--Offline','');
GenNum = find(strcmp(selection,names));
Plant.Generator(GenNum).Enabled = 1;
listbox_MakeList(hObject, eventdata, handles)

function pushbuttonEditData_Callback(hObject, eventdata, handles)
global figHandle
SetupBuilding
waitfor(figHandle)

function OptimOptions_Callback(hObject, eventdata, handles)
OptimOptions %all the options are loaded into the structure Plant.optimoptions

function pushbuttonRealTime_Callback(hObject, eventdata, handles)
global RealTime Virtual
Virtual = 0;
RealTime=1;
DispatchMPC(handles)

function pushbuttonDispatch_Callback(hObject, eventdata, handles)
global RealTime Virtual
Virtual = 1;
RealTime =0;
DispatchMPC(handles)


function DispatchMPC(handles)
global Plant Dispatch  DispatchWaitbar 
%% ---- %%%
RunHMPC
%%
waitfor(DispatchWaitbar)
Dispatch.Baseline = RunBaseline; %finish simulation by running baseline
closePorts
Plant.Dispatch = Dispatch;
set(handles.textDispatch,'string',num2str(round(sum(Plant.Dispatch.NetCost))))
if ~isfield(Plant.Dispatch,'Historical')
    days = max(1,round(Plant.Dispatch.Dispatch.Timestamp(end) - Plant.Dispatch.Dispatch.Timestamp(1)));
    set(handles.sliderDate2,'Min',1,'Max',2,'Value',1,'SliderStep',[1/days,1/days])
end
popupmenuGenerators_Callback(handles.popupmenuData,[],handles)
% ExportResult


%% re-do figure functions
function sliderDate1_Callback(hObject, eventdata, handles)
Plot2Axes(handles)
function sliderDate2_Callback(hObject, eventdata, handles)
popupmenuGenerators_Callback(handles.popupmenuData,eventdata,handles)
function checkboxEnergyStorageState_Callback(hObject, eventdata, handles)
popupmenuGenerators_Callback(hObject, eventdata, handles)
function popupmenuData_Callback(hObject, eventdata, handles)
Plot2Axes(handles)
function checkboxForecast_Callback(hObject, eventdata, handles)
Plot2Axes(handles)
%% Do nothing functions
function sliderZoom1_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
function sliderDate2_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function popupmenuPlant_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function sliderDate1_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
function popupmenuGenerators_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function popupmenuData_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function axes2_CreateFcn(hObject, eventdata, handles)
%% Menu functions
function uiAbout_Callback(hObject, eventdata, handles)
Message=['EAGERS Model. Developed at the Clean Energy Systems Integration Lab at Washington State University',...
    'to model and optimize micro-grid dispatch'];
msgbox(Message)
function uiProductHelp_Callback(hObject, eventdata, handles)
open('EAGERS User Guide.docx')
