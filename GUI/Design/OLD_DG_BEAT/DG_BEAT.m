function varargout = DG_BEAT(varargin)
% DG_BEAT MATLAB code for DG_BEAT.fig
%      DG_BEAT, by itself, creates a new DG_BEAT or raises the existing
%      singleton*.
%
%      H = DG_BEAT returns the handle to a new DG_BEAT or the handle to
%      the existing singleton*.
%
%      DG_BEAT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DG_BEAT.M with the given input arguments.
%
%      DG_BEAT('Property','Value',...) creates a new DG_BEAT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before DG_BEAT_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to DG_BEAT_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help DG_BEAT

% Last Modified by GUIDE v2.5 15-Nov-2015 15:46:23

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @DG_BEAT_OpeningFcn, ...
                   'gui_OutputFcn',  @DG_BEAT_OutputFcn, ...
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
% Added Comment

% --- Executes just before DG_BEAT is made visible.
function DG_BEAT_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to DG_BEAT (see VARARGIN)

% Choose default command line output for DG_BEAT
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes DG_BEAT wait for user response (see UIRESUME)
% uiwait(handles.figure1);

global Model_dir
set(gcf,'Name','DG-BEAT 2017.0.1')

%change to find project directory more reliably
projdir=fullfile(Model_dir, 'DesignProjects');
files=dir(fullfile(projdir,'*.mat'));
list=strrep({files.name},'.mat','');

global Project
if isempty(Project)
    set(handles.popupmenuProject,'string',list,'value',1)
else
    value=strmatch(Project.Name,list,'exact');
    set(handles.popupmenuProject,'string',list,'value',value)
end

plotlist={'Electric Demand'
        'Electric Demand Hist'
        'Electric without HVAC'
        'Electric Demand by Day'
        'Cooling Demand'
        'Cooling Hist'
        'Heat Demand'
        'Heat Demand Hist'};
set(handles.popupmenuAxes1,'string',plotlist,'value',1)
set(handles.popupmenuAxes2,'string',plotlist,'value',4)

set(handles.sliderZoom1,'Min',1,'Max',4,'Value',1,'SliderStep',[1/3,1/3])

popupmenuProject_Callback(handles.popupmenuProject,eventdata,handles,0)
sliderZoom1_Callback(handles.sliderZoom1, eventdata, handles)


    

% --- Outputs from this function are returned to the command line.
function varargout = DG_BEAT_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in popupmenuAxes1.
function popupmenuAxes_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuAxes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
tag=get(hObject,'tag');
axesnumstr=tag(end);
val=get(hObject,'value');
str=get(hObject,'string');
currentstr=str{val};

global Project Model_dir
Ts = 8760/length(Project.Building(1).DemandE);
A = gcbf;
if isempty(A)
    axes(findobj(gcf,'tag',['axes' axesnumstr]))
else axes(findobj(gcbf,'tag',['axes' axesnumstr]))
end
cla reset
set(gca,'tag',['axes' axesnumstr])
hold all

%Show result of test
name = char(Project.Building.Name);
%real building names

compdir=fullfile(Model_dir,'System Library', 'Buildings','RealBuildingData');
files=dir(fullfile(compdir,'*.mat'));
RealBuildNames=strrep({files.name},'.mat','');
if strncmp(name,'Multiple',8) %Refers to an aggragate of multiple building types
    set(handles.textBuild,'string','Multiple')
    set(handles.textClimate,'string','Multiple')
    set(handles.textVintage,'string','Multiple')
elseif max(strcmp(RealBuildNames,name))>0 %refers to a real building
    set(handles.textBuild,'string','Real Building')
    set(handles.textClimate,'string','Real Building')
    set(handles.textVintage,'string','Real Building')
else
    ind =strfind(name,'_');
    build = name(1:ind(1)-1);
    clim = name(ind(1):ind(2));
    vin = name(ind(2)+1:end);
    BuildType = {'Restaurant: full-service (sit down)';'Restaurant: quick-service (fast food)';'School: primary school';'School: secondary school';'Office: large office';'Office: medium office';'Office: small office';'Mid-rise apartment building';'Hospitality: large hotel';'Hospitality: small hotel/motel';'Health care: large hospital';'Health care: outpatient facility';'Retail: big-box, standalone retail store';'Retail: retail store located in a strip mall';'Retail: supermarket';'Unrefrigerated warehouse';};
    Climate = {'Miami (ASHRAE 1A)';'Houston (ASHRAE 2A)';'Phoenix (ASHRAE 2B)';'Atlanta (ASHRAE 3A)';'Las Vegas (ASHRAE 3B-Inland)';'Los Angeles (ASHRAE 3B-Coast)';'San Francisco (ASHRAE 3C)';'Baltimore (ASHRAE 4A)';'Albuquerque (ASHRAE 4B)';'Seattle (ASHRAE 4C)';'Chicago (ASHRAE 5A)';'Boulder (ASHRAE 5B)';'Minneapolis (ASHRAE 6A)';'Helena, MT (ASHRAE 6B)';'Duluth, MN (ASHRAE 7)';'Fairbanks, AK (ASHRAE 8)';};
    Vintage = {'2010 construction (ASHRAE 90.1-2010)';'2007 construction (ASHRAE 90.1-2007)';'2004 construction 90.1-2004';'“Post-1980” construction (ASHRAE 90.1-1989)';'“Pre-1980” construction';};
    vintageName = {'New2010';'New2007';'New2004';'Post1980';'Pre1980';};
    climateName = {'_1A_'; '_2A_'; '_2B_'; '_3A_'; '_3B_'; '_3B-Coast_'; '_3C_'; '_4A_'; '_4B_'; '_4C_'; '_5A_'; '_5B_'; '_6A_'; '_6B_'; '_7_'; '_8_';};
    buildTypeName = {'SDRest'; 'FFRest'; 'Sch-pri'; 'Sch-sec'; 'LgOff'; 'MdOff'; 'SmOff'; 'MRapt'; 'LgHotel'; 'SmHotel'; 'Hospital'; 'OutP'; 'Retail'; 'StMall'; 'SMarket'; 'ware';};
    set(handles.textBuild,'string',BuildType(find(ismember(buildTypeName, build))))
    set(handles.textClimate,'string',Climate(find(ismember(climateName, clim))))
    set(handles.textVintage,'string',Vintage(find(ismember(vintageName, vin))))
%     h=handles.axesImage;
%     cla(h)
%     Pic_dir=fullfile(Model_dir, 'graphics');
%     pFileName = strcat(Pic_dir,filesep ,char(buildTypeName(find(ismember(buildTypeName, build)))),'.jpg');
%     image(imread(pFileName),'Parent',h)
%     axis(h,'image')
%     axis(h,'off')
end

result = Project.Result.eOut;
cost = Project.Result.costOut;

costNum = [mean(cost.BaselineGridBill+cost.BaselineFuel);mean(cost.CostPerYear);mean(cost.FuelCosts);mean(cost.NewGridBill);];
if min(costNum)<1000
    costNum = round(costNum);
elseif min(costNum)<1e4
    costNum = 10*round(costNum/10);
elseif min(costNum)<1e5
    costNum = 100*round(costNum/100);
elseif min(costNum)<1e6
    costNum = 1e3*round(costNum/1e3);
elseif min(costNum)<1e7
    costNum = 1e4*round(costNum/1e4);
else costNum = 1e5*round(costNum/1e5);
end
ResultsLabel = {strcat('CHP size (kW):');
            strcat('CHP capacity factor:');
            strcat('Self Generation proportion:');
            strcat('Heating met by CHP (%):');
            strcat('Baseline Annual Cost ($):');
            strcat('With DG Annual Cost ($):');
            strcat('DG Fuel Costs ($): ');
            strcat('With DG Grid Cost ($):');};
        
Results = {strcat(num2str(round(result.SysSize)));
            strcat(num2str(result.Total_Electricity_produced_kWh/(result.SysSize*8760),2));
            strcat(num2str(result.SelfGen,2));
            strcat(num2str(round(100*result.Total_DG_Heat_out_kWh/(result.Total_DG_Heat_out_kWh+ result.Total_peak_burner_heat_out_kWh))));
            strcat(num2str(round(costNum(1))));
            strcat(num2str(round(costNum(2))));
            strcat(num2str(round(costNum(3))));
            strcat(num2str(round(costNum(4))));};
set(handles.textResultsLabel,'string',ResultsLabel)
set(handles.textResults,'string',Results)
BuildingSummaryLabel = {strcat('Max Load (kW):');
                        strcat('Min Load (kW):');
                        strcat('Ave Load (kW):');
                        strcat('% Electricity for Cooling:')};
BuildingSummary = {strcat(num2str(round(max(Project.Building.DemandE))));
                    strcat(num2str(round(min(Project.Building.DemandE))));
                    strcat(num2str(round(mean(Project.Building.DemandE))));
                    strcat(num2str(100*sum(Project.Building.CoolingElectricalLoad)/sum(Project.Building.DemandE),2));
                    };
set(handles.textBuildingStats,'string',BuildingSummary)
set(handles.textBuildingLabel,'string',BuildingSummaryLabel)

Zoom = max(1,round(get(handles.sliderZoom1,'Value')));
StartDate = max(1,floor(get(handles.sliderDate,'Value')));
month = [0 31 28 31 30 31 30 31 31 30 31 30 31];
monthDays = [0 31 59 90 120 151 181 212 243 273 304 334 365];
monthLabel = ['January  '; 'February '; 'March    '; 'April    '; 'May      '; 'June     '; 'July     '; 'August   '; 'September'; 'October  ';'November ' ;'December ';];
if Zoom == 1;
    day1 = 1;
    lastDay = 365;
    plotAxis = datenum([2014*ones(12,1) (1:12)' 1*ones(12,1) 0*ones(12,1) 0*ones(12,1) 0*ones(12,1)]);
    xlabel('Date')
elseif Zoom == 2;
    day1 = 1+sum(month(1:StartDate));
    lastDay = sum(month(1:StartDate+1));
    days = month(StartDate+1);
    plotAxis = datenum([2014*ones(days,1) StartDate*ones(days,1) (1:days)' 0*ones(days,1) 0*ones(days,1) 0*ones(days,1)]);
	xlabel(monthLabel(StartDate,:))
elseif Zoom == 3;
    day1 = 1+7*(StartDate-1);
    lastDay = 1+7*StartDate;
    StartMonth = find(monthDays>day1,1) - 1;
    StartDay = day1-monthDays(StartMonth);
    plotAxis = datenum([2014*ones(8,1) StartMonth*ones(8,1)  (StartDay:StartDay+7)'  0*ones(8,1)  0*ones(8,1) 0*ones(8,1)]);
    if StartDay+6<=month(StartMonth+1)
        xlabel(monthLabel(StartMonth,:))
    else 
        xlabel(strcat(monthLabel(StartMonth,:), ' / ', monthLabel(StartMonth+1,:)))
    end
elseif Zoom == 4
    day1 = StartDate;
    lastDay = StartDate;
    StartMonth = find(monthDays>day1,1) - 1;
    StartDay = day1-monthDays(StartMonth);
    plotAxis = datenum([2014*ones(24,1)  StartMonth*ones(24,1)  StartDay*ones(24,1) (1:24)' 0*ones(24,1) 0*ones(24,1)]);
    xlabel(strcat(['Hours of ', monthLabel(StartMonth,:), num2str(StartDay)]))
end
X = datenum(2014,1,1,0,0,0)+((day1-1)+(Ts/24):(Ts/24):lastDay);
switch currentstr
    case 'Heat Demand'
        PlotType = 'plot';
        Y = Project.Building(1).DemandH(1+24/Ts*(day1-1):24/Ts*lastDay);
    case 'Heat Demand Hist'
        PlotType = 'hist';
        Y = Project.Building(1).DemandH(1+24/Ts*(day1-1):24/Ts*lastDay);
    case 'Electric Demand'
        PlotType = 'plot';
        Y = Project.Building(1).DemandE(1+24/Ts*(day1-1):24/Ts*lastDay);
    case 'Electric Demand by Day'
        PlotType = 'pcolor';
        xgrid=floor(X(1)):floor(X(end));
        ygrid=Ts:Ts:24;
        z=NaN*ones(length(xgrid),length(ygrid));
        for i=1:lastDay-day1+1
            z(i,:)=Project.Building(1).DemandE(1+24/Ts*(i-1+day1-1):(i+day1-1)*24/Ts);
        end
    case 'Electric Demand Hist'
        PlotType = 'hist';
        Y = Project.Building(1).DemandE(1+24/Ts*(day1-1):24/Ts*lastDay);
    case 'Cooling Demand'
        PlotType = 'plot';
        Y = Project.Building(1).DemandC(1+24/Ts*(day1-1):24/Ts*lastDay);
    case 'Cooling Hist'
        PlotType = 'hist';
        Y = Project.Building(1).DemandC(1+24/Ts*(day1-1):24/Ts*lastDay);
    case 'Electric without HVAC'
        PlotType = 'plot';
        Y = Project.Building(1).DemandE(1+24/Ts*(day1-1):24/Ts*lastDay)-Project.Building(1).CoolingElectricalLoad(1+24/Ts*(day1-1):24/Ts*lastDay);
end
if strcmp(PlotType,'plot')
    plot(X,Y);
    ylabel('Demand (kW)') 
elseif strcmp(PlotType,'hist')
    hist(Y);
    xlabel('Demand (kW)')
    ylabel('Count')
elseif strcmp(PlotType,'pcolor')
    h=pcolor(xgrid,ygrid,z');
    set(h,'edgecolor','none')
    ylabel('hour')
    colorbar
    h = colorbar;
    ylabel(h,'(kW)');
    ylim([1 24]) 
elseif strcmp(PlotType,'surface')
    surface(S,Y,Z,'facecolor','interp','edgecolor','none')
    xlabel('FC Size [kW]')
    ylabel('Production Quantity')
    zlabel('Cost per kW [$]')
    colorbar
    view(206,12)
end
if strcmp(PlotType,'plot') || strcmp(PlotType,'pcolor')
    set(gca,'xtick',plotAxis) 
    if Zoom ==1
        datetick('x','mmmdd','keepticks')
    elseif Zoom ==2
        datetick('x','dd','keepticks')
    elseif Zoom ==3
        datetick('x','dd','keepticks')
    elseif Zoom ==4
        datetick('x','HH','keepticks')
    end
end
str = get(hObject,'string');
elec1 = strcmp(str(get(handles.popupmenuAxes1,'value')),'Electric Demand') || strcmp(str(get(handles.popupmenuAxes1,'value')),'Electric without HVAC');
elec2 = strcmp(str(get(handles.popupmenuAxes2,'value')),'Electric Demand') || strcmp(str(get(handles.popupmenuAxes2,'value')),'Electric without HVAC');
if elec1 || elec2
    set(handles.checkboxGeneration,'enable','on')
else set(handles.checkboxGeneration,'enable','off')
end
if isfield(Project.System,'TES')
    set(handles.checkboxTESshift,'enable','on')
else set(handles.checkboxTESshift,'enable','off')
end
if isfield(Project.System,'Battery')
    set(handles.checkboxBattery,'enable','on')
else set(handles.checkboxBattery,'enable','off')
end
if isfield(Project.System,'PV') || isfield(Project.System,'Wind')
    set(handles.checkboxRenewable,'enable','on')
else set(handles.checkboxRenewable,'enable','off')
end
if strcmp(str(get(hObject,'value')),'Electric Demand') || strcmp(str(get(hObject,'value')),'Electric without HVAC')
    if get(handles.checkboxGeneration,'Value')==1
        genPow = result.GenPower(1+24/Ts*(day1-1):24/Ts*lastDay);
        hold on
        plot(X,genPow,'r');
        hold off
    end
    if get(handles.checkboxTESshift,'Value')==1 && isfield(Project.System,'TES')
        hold on
        plot(X,result.DemandEshift(1+24/Ts*(day1-1):24/Ts*lastDay),'g');
        hold off
    end
    if get(handles.checkboxBattery,'Value')==1 && isfield(Project.System,'Battery')
        hold on
        plot(X,result.DemandEshift2(1+24/Ts*(day1-1):24/Ts*lastDay),'m');
        hold off
    end
end


% --- Executes during object creation, after setting all properties.
function popupmenuAxes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuAxes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function popupmenuAxes2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuAxes2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in popupmenuProject.
function popupmenuProject_Callback(hObject, eventdata, handles, skipreload)
% hObject    handle to popupmenuProject (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Project Model_dir
ProjectList=get(hObject,'string');
ProjectVal=get(hObject,'value');
if nargin ==3 || skipreload == 0
    projdir=fullfile(Model_dir,'DesignProjects');
    if ~isempty(ProjectVal)
        projName=fullfile(projdir,ProjectList{ProjectVal});
    else
        projName=fullfile(projdir,ProjectList{1});
    end
    load (projName)
else
    Project.Result = runAnalyses(Project);
end

if ~isfield(Project,'State')
    Project.State = 'California';
end
if isfield(Project.Utilities.NatGas,'Type')
    load(fullfile(Model_dir,'System Library','NatGas','Natural Gas Prices Default'))
    Project.Utilities.NatGas = component;
end
%build list for Setup Listbox
fields=fieldnames(Project);
str={};

for i=1:length(fields)
    if isstruct(Project.(fields{i})) && (1-strcmp(fields{i},'Result'))
          str{end+1}=[fields{i} ':' Project.(fields{i}).Name];
          fields2=fieldnames(Project.(fields{i}));
          if strcmp(fields{i},'Building') && length(Project.(fields{i}))>1
              if ~strcmp(str{end},'Building:Multiple')
                str{end+1}='Building:Multiple';
              end
          else
              for j=1:length(fields2)
                  if isstruct(Project.(fields{i}).(fields2{j}))
                      %only use structures that have Name fields otherwise
                      %structure is part of the component data and not the top level for that component

                      if isfield(Project.(fields{i}).(fields2{j}),'Name') 
                          for k=1:length(Project.(fields{i}).(fields2{j}))
                             str{end+1}=['   '   fields{i} '.' fields2{j} ':' Project.(fields{i}).(fields2{j})(k).Name];
                          end
                      end
                  end
              end
          end
    end
end
set(handles.listboxSetup,'string',str,'value',1)
listboxSetup_Callback(handles.listboxSetup,eventdata,handles)
popupmenuAxes_Callback(handles.popupmenuAxes1,eventdata,handles)
popupmenuAxes_Callback(handles.popupmenuAxes2,eventdata,handles)

% --- Executes during object creation, after setting all properties.
function popupmenuProject_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuProject (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pushbuttonProject.
function pushbuttonProject_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonProject (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Project Model_dir
name = char(Project.Name);
projdir=fullfile(Model_dir,'DesignProjects');

[f,p]=uiputfile('*.mat','Save Project As...',fullfile(projdir, name));
if f==0; return; end
saveAs(Project,f,p)

DG_BEAT_OpeningFcn(handles.pushbuttonProject,eventdata,handles)

function saveAs(Project,filename,saveDirectory)
Project.id=str2double(datestr(now,'yymmddHHMMSS'));
Project.Name=strrep(filename,'.mat','');
save([saveDirectory,filename],'Project')

% --- Executes on selection change in listboxSetup.
function listboxSetup_Callback(hObject, eventdata, handles)
% hObject    handle to listboxSetup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%get the component selected
list=get(handles.listboxSetup,'string');
val=get(handles.listboxSetup,'value');
compSel=getCurrentComponent(list,val);
inddot=strfind(compSel,'.');
if isempty(inddot)
    inddot = length(compSel)+1;
end
oncomponents={'System'};
oncomponents2={'Renewable'};
if strcmp(compSel(1:inddot(1)-1),oncomponents)
   set(handles.pushbuttonAdd,'enable','on')
   set(handles.pushbuttonRemove,'enable','on')
elseif strcmp(compSel(1:inddot(1)-1),oncomponents2)
   set(handles.pushbuttonAdd,'enable','on')
   set(handles.pushbuttonRemove,'enable','on')
else
    set(handles.pushbuttonAdd,'enable','off')
    set(handles.pushbuttonRemove,'enable','off')
end  

% --- Executes during object creation, after setting all properties.
function listboxSetup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listboxSetup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over listboxSetup.
function listboxSetup_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to listboxSetup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global SYSINDEX Project
listnum=get(hObject,'value');
str=get(hObject,'string');
selection=str{listnum};
if ~isempty(strfind(selection,'NatGas'))
    NatGasSetup()
elseif ~isempty(strfind(selection,'Economic'))
    Economics()
elseif ~isempty(strfind(selection,'System.CHP'))
    ind=find(strcmp(validatestring('System.CHP',strtrim(str)),strtrim(str))==1);
    SYSINDEX=find(ind==listnum);
    FuelCellSetup()
elseif ~isempty(strfind(selection,'System.Chiller')) 
    ind=strmatch('System.Chiller',strtrim(str));
    SYSINDEX=find(ind==listnum);
    ChillerSetup()
elseif ~isempty(strfind(selection,'System.TES'))
    ind=strmatch('System.TES',strtrim(str));
    SYSINDEX=find(ind==listnum);
    ThermalStorageSetup()
elseif ~isempty(strfind(selection,'System.Battery'))
    ind=strmatch('System.Battery',strtrim(str));
    SYSINDEX=find(ind==listnum);
    BatterySetup()
elseif ~isempty(strfind(selection,'Grid'))
    ModifyElectricGridRate()   
elseif ~isempty(strfind(selection,'Renewable.Solar'))
    ind=strmatch('Renewable.Solar',strtrim(str));
    SYSINDEX=find(ind==listnum);
    SolarSetup()
elseif ~isempty(strfind(selection,'Renewable.Wind'))
    ind=strmatch('Renewable.Wind',strtrim(str));
    SYSINDEX=find(ind==listnum);
    WindSetup()
else
	uiresume
end
uiwait
Project.Result = runAnalyses(Project);
popupmenuAxes_Callback(handles.popupmenuAxes1,eventdata,handles)
popupmenuAxes_Callback(handles.popupmenuAxes2,eventdata,handles)


% --- Executes on button press in pushbuttonAdd.
function pushbuttonAdd_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonAdd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
h=ComponentSelector();
uiwait(h)
popupmenuProject_Callback(handles.popupmenuProject,eventdata,handles,1);


% --- Executes on button press in pushbuttonChange.
function pushbuttonChange_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonChange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Project Model_dir
list = get(handles.listboxSetup,'string');
change =  get(handles.listboxSetup,'value');
componentSelected=getCurrentComponent(list,change);
compindex = 0;
for i = 1:1:change-1
    if strcmp(getCurrentComponent(list,i),componentSelected)
        compindex = compindex+1;
    end
end
componentSel=strrep(componentSelected,'System.','');
componentSel=strrep(componentSel,'Utilities.','');
componentSel=strrep(componentSel,'Renewable.','');
componentSel=strrep(componentSel,'Utilities','Grid');
componentSel=strrep(componentSel,'TES','ThermalStorage');
if  strcmp(componentSel, 'System')
    componentSel = 'CHP';
    compindex = 0;
end
if strcmp(componentSel, 'Building')
    buildingChooser();
    uiwait(gcf)
elseif strcmp(componentSel,'NatGas')
    NatGasSetup()
    uiwait
else
    compdir=fullfile(Model_dir,'System Library', componentSel);
    files=dir(fullfile(compdir,'*.mat'));
    list=strrep({files.name},'.mat','');
    [selection, ok]=listdlg('ListString',list);
    if ok
        selection=selection(1);
        componentSel=strrep(componentSel,'ThermalStorage','TES');
        load(strcat(compdir,filesep,list{selection})) %load component structure
        if strncmp(componentSelected,'System',min(length(componentSelected),6))
            Project.System.(componentSel)(compindex+1)=component;
        elseif strncmp(componentSelected,'Utilities',min(length(componentSelected),9))
            Project.Utilities.(componentSel) = component;
        elseif strncmp(componentSelected,'Renewable',min(length(componentSelected),9))
            Project.Renewable.(componentSel) = component;
        else Project.(componentSel) = component;
        end
    end
end
popupmenuProject_Callback(handles.popupmenuProject, eventdata, handles,1)


% --- Executes on button press in pushbuttonRemove.
function pushbuttonRemove_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonRemove (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
list=get(handles.listboxSetup,'string');
val=get(handles.listboxSetup,'value');
componentSelected=getCurrentComponent(list,val);
global Project
answer=questdlg({'Are you sure you want to remove:'
    componentSelected});
if strcmp(answer,'Yes')
    %only remove last field
    ind=strfind(componentSelected,'.');
    if isempty(ind) && strcmp(componentSelected,'Renewable')
        Project = rmfield(Project,'Renewable');
    else
        str=componentSelected;
        category = str(1:ind(1)-1);
        currentfield=str(ind(1)+1:end);
        ind2struct=1;
        prev = 1;
        while ~isempty(strfind(getCurrentComponent(list,val-prev),'.'))% && (val-ind2struct>0)
            if strcmp(str,getCurrentComponent(list,val-prev));
                ind2struct=ind2struct+1;
            end
            prev = prev+1;
        end
        j=0;
        for i = 1:1:length(Project.(category).(currentfield))
            if i~=ind2struct
                j = j+1;
                COMPONENT(j) = Project.(category).(currentfield)(i);
            end
        end
        Project.(category)= rmfield(Project.(category),currentfield);
        if j>0
            Project.(category).(currentfield) = COMPONENT;
        end
        if isfield(Project,'Renewable')
            if length(fieldnames(Project.Renewable))==1
                Project = rmfield(Project,'Renewable');
            end
        end
    end
    skipreload=1;
    popupmenuProject_Callback(handles.popupmenuProject,eventdata,handles,skipreload);
end


% --- Executes on slider movement.
function sliderZoom1_Callback(hObject, eventdata, handles)
% hObject    handle to sliderZoom1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Zoom = round(get(hObject,'Value'));
if Zoom == 1
    steps = 2;
elseif Zoom == 2
    steps = 12;
elseif Zoom == 3
    steps = 52;
elseif Zoom == 4
    steps = 365;
end
set(handles.sliderDate,'Min',1)
set(handles.sliderDate,'Max',steps)
set(handles.sliderDate,'Value',1,'SliderStep',[1/(steps-1),1/(steps-1)])
popupmenuAxes_Callback(handles.popupmenuAxes1,eventdata,handles)
popupmenuAxes_Callback(handles.popupmenuAxes2,eventdata,handles)

% --- Executes during object creation, after setting all properties.
function sliderZoom1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderZoom1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function sliderDate_Callback(hObject, eventdata, handles)
% hObject    handle to sliderDate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
popupmenuAxes_Callback(handles.popupmenuAxes1,eventdata,handles)
popupmenuAxes_Callback(handles.popupmenuAxes2,eventdata,handles)

% --- Executes during object creation, after setting all properties.
function sliderDate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderDate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes on button press in checkboxGeneration.
function checkboxGeneration_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxGeneration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
popupmenuAxes_Callback(handles.popupmenuAxes1,eventdata,handles)
popupmenuAxes_Callback(handles.popupmenuAxes2,eventdata,handles)


% --- Executes on button press in checkboxTESshift.
function checkboxTESshift_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxTESshift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
popupmenuAxes_Callback(handles.popupmenuAxes1,eventdata,handles)
popupmenuAxes_Callback(handles.popupmenuAxes2,eventdata,handles)

% --- Executes on button press in checkboxBattery.
function checkboxBattery_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxBattery (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
popupmenuAxes_Callback(handles.popupmenuAxes1,eventdata,handles)
popupmenuAxes_Callback(handles.popupmenuAxes2,eventdata,handles)

% --- Executes on button press in checkboxRenewable.
function checkboxRenewable_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxRenewable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
popupmenuAxes_Callback(handles.popupmenuAxes1,eventdata,handles)
popupmenuAxes_Callback(handles.popupmenuAxes2,eventdata,handles)

% --- Executes on button press in pushbuttonResults.
function pushbuttonResults_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonResults (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
currentProject=get(handles.popupmenuProject, 'value');
ViewResults2(currentProject)


% --- Executes on button press in pushbuttonOptSize.
function pushbuttonOptSize_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonOptSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ScaleSystemComponents(100,100,100,120,1,2,1,2,2)
popupmenuProject_Callback(handles.popupmenuProject,eventdata,handles,1);

% --- Executes on button press in pushbuttonFCcostSize.
function pushbuttonFCcostSize_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonFCcostSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Project
if isfield(Project.System,'TES')
    autoSizing('cost','TES')
end
autoSizing('cost','CHP')
if isfield(Project.System,'Battery')
    autoSizing('cost','Battery',3)
end
popupmenuProject_Callback(handles.popupmenuProject,eventdata,handles,1);

% --- Executes on button press in pushbuttonEmSize.
function pushbuttonEmSize_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonEmSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
list={'CO2'; 'SO2'; 'NOx';};
    [selection, ok]=listdlg('ListString',list);
    if ok
        selection=list{selection(1)};
        switch selection
            case 'CO2'
                autoSizing('CO2','CHP')
            case 'SO2'
                autoSizing('SO2','CHP')
            case 'NOx'
                autoSizing('NOx','CHP')
        end
    end
popupmenuProject_Callback(handles.popupmenuProject,eventdata,handles,1);

function componentSelected=getCurrentComponent(list,val)
componentSelected=strtrim(strtok(list{val},':'));

% --------------------------------------------------------------------
function uiAbout_Callback(hObject, eventdata, handles)
% hObject    handle to uiAbout (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Message=['Distributed Generation Build-Out and Economic assesment Tool (DG-BEAT)' ,...
    ' Designed to optimize the use of distributed energy technologies at commercial buildings.'];
msgbox(Message)

% --------------------------------------------------------------------
function uiProductHelp_Callback(hObject, eventdata, handles)
% hObject    handle to uiProductHelp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
open('Stationary Fuel Cell Model User Guide.docx')

% --------------------------------------------------------------------
function uIresults_Callback(hObject, eventdata, handles)
% hObject    handle to uIresults (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
currentProject=get(handles.popupmenuProject, 'value');
% ViewResults(currentProject)
ViewResults2(currentProject)

% --------------------------------------------------------------------
function uImultiBuildResult_Callback(hObject, eventdata, handles)
% hObject    handle to uImultiBuildResult (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
MultiBuildResult()

% --------------------------------------------------------------------
function uIsensitivityResults_Callback(hObject, eventdata, handles)
% hObject    handle to uIsensitivityResults (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
SensitivityResult()

% --------------------------------------------------------------------
function uImultiBuild_Callback(hObject, eventdata, handles)
% hObject    handle to uImultiBuild (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
MultiBuildingAnalysis()

% --------------------------------------------------------------------
function uIsensitivity_Callback(hObject, eventdata, handles)
% hObject    handle to uIsensitivity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Sensitivity()

% --------------------------------------------------------------------
function uIrealBuild_Callback(hObject, eventdata, handles)
% hObject    handle to uIrealBuild (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Validation()

% --------------------------------------------------------------------
function uIchangeBuild_Callback(hObject, eventdata, handles)
% hObject    handle to uIchangeBuild (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
buildingChooser();
uiwait
popupmenuProject_Callback(handles.popupmenuProject, eventdata, handles,1)

% --------------------------------------------------------------------
function uIcontrol_Callback(hObject, eventdata, handles)
% hObject    handle to uIcontrol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Project Model_dir
Control_dir=fullfile(Model_dir, 'System Library' , 'Control'); 
files=dir(fullfile(Control_dir,'*.mat'));
list=strrep({files.name},'.mat','');
[selection, ok]=listdlg('ListString',list);
if ok
    load(list{selection(1)}) %load component structure
    Project.Control = component;
end
popupmenuProject_Callback(handles.popupmenuProject, eventdata, handles,1)

% --------------------------------------------------------------------
function uIeconomics_Callback(hObject, eventdata, handles)
% hObject    handle to uIeconomics (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Economics()

% --------------------------------------------------------------------
function uIgrid_Callback(hObject, eventdata, handles)
% hObject    handle to uIgrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ModifyElectricGridRate() 

% --------------------------------------------------------------------
function uiNationalSurvey_Callback(hObject, eventdata, handles)
% hObject    handle to uiNationalSurvey (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
NationalSurvey()

% --------------------------------------------------------------------
function uiNationalSurveyResults_Callback(hObject, eventdata, handles)
% hObject    handle to uiNationalSurvey (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
NationalSurveyResults()

% --------------------------------------------------------------------
function FleetAnalysis_Callback(hObject, eventdata, handles)
% hObject    handle to FleetAnalysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
FleetAnalysis()

% --------------------------------------------------------------------
function FleetResult_Callback(hObject, eventdata, handles)
% hObject    handle to FleetResult (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
FleetResult()


% --- Executes on button press in pushbuttonHelpMe.
function pushbuttonHelpMe_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonHelpMe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
HelpScreen1()
