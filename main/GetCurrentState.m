function [State,CoGen,OnOff,Input] = GetCurrentState
global MeasurePort MeasurePortIndex

for i = 1:1:length(MeasurePort)
    fopen(MeasurePort(i));%open all ports
end
%initialize variables
State = zeros(1,length(MeasurePortIndex));
CoGen = zeros(1,length(MeasurePortIndex));
OnOff = zeros(1,length(MeasurePortIndex));
Input = zeros(1,length(MeasurePortIndex));
StateJ =[];
StateI =[];
CoGenJ =[];
CoGenI =[];
OnOffJ =[];
OnOffI =[];
InputJ =[];
InputI =[];
%organize what each datum goes to
for i = 1:1:length(MeasurePortIndex)
    if MeasurePortIndex(i,3)~=0
        StateJ(end+1) = MeasurePortIndex(i,3);
        StateI(end+1) = i;
        if MeasurePortIndex(i,4)~=0
            CoGenJ(end+1) = MeasurePortIndex(i,4);
            CoGenI(end+1) = i;
        end
    else
        if MeasurePortIndex(i,4)~=0
            StateJ(end+1) = MeasurePortIndex(i,4);
            StateI(end+1) = i;
        end
    end
    if MeasurePortIndex(i,2)~=0
        InputJ(end+1) = MeasurePortIndex(i,2);
        InputI(end+1) = i;
    end
    if MeasurePortIndex(i,1)~=0
        OnOffJ(end+1) = MeasurePortIndex(i,1);
        OnOffI(end+1) = i;
    end
    %% battery must be measured with 0 = SOC at maximum discharge and max value is battery size* max discharge (only the usable kWh of the battery)
end
Data = zeros(length(MeasurePort),1);
for i = 1:1:length(MeasurePort)
    Data(i) = round(str2double(char(fread(MeasurePort(i))))); %read all data
end
for i = 1:1:length(MeasurePort) 
    fclose(MeasurePort(i)); %close all ports
end
%organize data into different vectors
State(StateI) = Data(StateJ);
CoGen(CoGenI) = Data(CoGenJ);
OnOff(OnOffI) = Data(OnOffJ);
Input(InputI) = Data(InputJ);
%% --- Logs Data during a real Test
LogData(State,CoGen,OnOff,Input)
%%%