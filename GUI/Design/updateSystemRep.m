function updateSystemRep(hObject, eventdata, handles)
%% System representation (pushbuttons need to appear/disappear if they exist)
global Plant

% Hide everything except grid and AC/DC conversion
set(handles.pushbuttonHeaterInSys,'Visible','off')
set(handles.textAirHeater_H1,'Visible','off')
set(handles.pushbuttonTES1,'Visible','off')
set(handles.textTES1_H1,'Visible','off')
set(handles.pushbuttonHeatingDemands,'Visible','off')
set(handles.pushbuttonICE_mGT,'Visible','off')
set(handles.textICE_mGT_E,'Visible','off')
set(handles.textICE_mGT_H2,'Visible','off')
set(handles.textICE_mGT_H3,'Visible','off')
set(handles.pushbuttonSolarThermalInSys,'Visible','off')
set(handles.textSolarThermal_H2,'Visible','off')
set(handles.pushbuttonWaterHeaterInSys,'Visible','off')
set(handles.textWaterHeater_H2,'Visible','off')
set(handles.pushbuttonTES2,'Visible','off')
set(handles.textTES2_H2,'Visible','off')
set(handles.pushbuttonFuelCell,'Visible','off')
set(handles.textFuelCell_E,'Visible','off')
set(handles.textFuelCell_H1,'Visible','off')
set(handles.textFuelCell_H2,'Visible','off')
set(handles.pushbuttonHotWaterDemands,'Visible','off')
set(handles.pushbuttonAbChillerInSys,'Visible','off')
set(handles.textAbChill_E,'Visible','off')
set(handles.textAbChill_C1,'Visible','off')
set(handles.textAbChill_C2,'Visible','off')
set(handles.pushbuttonChillerInSys,'Visible','off')
set(handles.textChill_E,'Visible','off')
set(handles.textChill_C1,'Visible','off')
set(handles.textChill_C2,'Visible','off')
set(handles.pushbuttonTES3,'Visible','off')
set(handles.textTES3_C1,'Visible','off')
set(handles.pushbuttonBatteryInSys,'Visible','off')
set(handles.textBattery_E,'Visible','off')
set(handles.pushbuttonCoolingDemands,'Visible','off')
set(handles.pushbuttonWindInSys,'Visible','off')
set(handles.textWind_E,'Visible','off')
set(handles.pushbuttonSolarSterlingInSys,'Visible','off')
set(handles.pushbuttonSolarPVInSys,'Visible','off')
set(handles.textSolarPV_E,'Visible','off')
% Show AC/DC conversion by default
set(handles.pushbuttonACDC,'Visible','on')
set(handles.textACDC_E,'Visible','on')

% Selectively show existing component categories
% Does not handle:
%   Absorption Chiller
%   Solar Stirling
%   Solar Thermal
%   TES 1
%   Wind
nG = length(Plant.Generator);
for i = 1:1:nG
    sys = Plant.Generator(i);
    if strcmp(sys.Type,'Utility') && strcmp(sys.Source,'Electricity')
        set(handles.pushbuttonGrid,'Visible','on')
        set(handles.textGrid,'Visible','on')
    elseif strcmp(sys.Type,'Heater')
        set(handles.pushbuttonHeaterInSys,'Visible','on')
        set(handles.textAirHeater_H1,'Visible','on')
    elseif strcmp(sys.Type,'Thermal Storage')
        if strcmp(sys.Source,'Heat')
            set(handles.pushbuttonTES2,'Visible','on')
            set(handles.textTES2_H2,'Visible','on')
        elseif strcmp(sys.Source,'Cooling')
            set(handles.pushbuttonTES3,'Visible','on')
            set(handles.textTES3_C1,'Visible','on')
        end
    elseif strcmp(sys.Type,'CHP Generator')
        if sys.VariableStruct.isFuelCell
            set(handles.pushbuttonFuelCell,'Visible','on')
            set(handles.textFuelCell_E,'Visible','on')
            set(handles.textFuelCell_H1,'Visible','on')
            set(handles.textFuelCell_H2,'Visible','on')
        else
            set(handles.pushbuttonICE_mGT,'Visible','on')
            set(handles.textICE_mGT_E,'Visible','on')
            set(handles.textICE_mGT_H2,'Visible','on')
        end
    elseif strcmp(sys.Type,'Solar')
        set(handles.pushbuttonSolarPVInSys,'Visible','on')
        set(handles.textSolarPV_E,'Visible','on')
    elseif strcmp(sys.Type,'Boiler')
        set(handles.pushbuttonWaterHeaterInSys,'Visible','on')
        set(handles.textWaterHeater_H2,'Visible','on')
    elseif strcmp(sys.Type,'Chiller')
        set(handles.pushbuttonChillerInSys,'Visible','on')
        set(handles.textChill_E,'Visible','on')
        set(handles.textChill_C1,'Visible','on')
        set(handles.textChill_C2,'Visible','on')
    elseif strcmp(sys.Type,'Electric Storage')
        set(handles.pushbuttonBatteryInSys,'Visible','on')
        set(handles.textBattery_E,'Visible','on')
    elseif strcmp(sys.Type,'Electric Generator')
        set(handles.pushbuttonICE_mGT,'Visible','on')
        set(handles.textICE_mGT_E,'Visible','on')
        set(handles.textICE_mGT_H2,'Visible','on')
    end
end

% Check whether hot/cold demands should be shown
if strcmp(get(handles.textICE_mGT_H2,'Visible'),'on') || ...
        strcmp(get(handles.textSolarThermal_H2,'Visible'),'on') || ...
        strcmp(get(handles.textWaterHeater_H2,'Visible'),'on') || ...
        strcmp(get(handles.textTES2_H2,'Visible'),'on') || ...
        strcmp(get(handles.textFuelCell_H2,'Visible'),'on')
    set(handles.pushbuttonHotWaterDemands,'Visible','on')
end
if strcmp(get(handles.textAirHeater_H1,'Visible'),'on') || ...
        strcmp(get(handles.textTES1_H1,'Visible'),'on') || ...
        strcmp(get(handles.textFuelCell_H1,'Visible'),'on')
    set(handles.pushbuttonHeatingDemands,'Visible','on')
end
if strcmp(get(handles.textAbChill_C1,'Visible'),'on') || ...
        strcmp(get(handles.textChill_C1,'Visible'),'on') || ...
        strcmp(get(handles.textTES3_C1,'Visible'),'on') || ...
        strcmp(get(handles.textTES2_H2,'Visible'),'on')
    set(handles.pushbuttonCoolingDemands,'Visible','on')
end