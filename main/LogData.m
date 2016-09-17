function LogData(State,CoGen,OnOff,Input)
%logs data specific to the plant being tested
global DataLog nLog FanPortRead AmbTempRead
nLog = nLog+1;
DataLog.Timestamp(nLog) = now;
%% %%% Read Fans
fopen(FanPortRead);
DataLog.FanTherm(nLog) =(str2double(char(fread(FanPortRead))));
fclose(FanPortRead);
%% Read Ambient Temperature
fopen(AmbTempRead);
DataLog.AmbTemp(nLog) =(str2double(char(fread(AmbTempRead))));
fclose(AmbTempRead);


DataLog.mGTpow(nLog) = State(5);
DataLog.mGTheat(nLog) = CoGen(5);
DataLog.mGTstate(nLog) = OnOff(5);
DataLog.mGTfuel(nLog) = Input(5);

DataLog.ICEpow(nLog) = State(4);
DataLog.ICEheat(nLog) = CoGen(4);
DataLog.ICEstate(nLog) = OnOff(4);
DataLog.ICEfuel(nLog) = Input(4);

DataLog.TES_SOC(nLog) = State(3); %kJ
%% verify reading is going well
disp(strcat('Current SOC measurement (MJ):',num2str(State(3)/1e3)))