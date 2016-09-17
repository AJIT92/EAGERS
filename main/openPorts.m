function openPorts
%% load/open comunication ports
global Plant CommPort CommPortIndex CommPortOnOff CommPortOnOffIndex MeasurePort MeasurePortIndex FanPortRead  AmbTempRead FanPortWrite
CommPort =[];
CommPortOnOff =[];
MeasurePort =[];
MeasurePortIndex = zeros(length(Generators),4);

if RealTime %open communication ports for turning generators on/off
    FanPortRead = udp('192.168.0.8', 31900, 'LocalPort',41900); %fan thermal power [kW]
    AmbTempRead = udp('192.168.0.8', 32000, 'LocalPort',42000); %Ambient Temperature [C]
    FanPortWrite = udp('192.168.0.8', 51500, 'LocalPort',61500); %set fan thermal demand (kW)
    fopen(FanPortWrite);

    Names = {'OnOff';'Input';'Electric';'Thermal';};
    for i= 1:1:length(Plant.Generator)
        if  isfield(Plant.Generator(i).VariableStruct,'Comm')
            a = Plant.Generator(i).VariableStruct.Comm.OnOff;
            b = Plant.Generator(i).VariableStruct.Comm.Set;
            if a~=0 
                port = udp('192.168.0.8', a, 'LocalPort',a+10000); 
                fopen(port);
                if isempty(CommPortOnOff)
                    CommPortOnOff = port;
                    CommPortOnOffIndex(i) = 1;
                else CommPortOnOff(end+1) = port;
                    CommPortOnOffIndex(i) = length(CommPortOnOff);
                end
            end
            if b~=0 
                port = udp('192.168.0.8', b, 'LocalPort',b+10000); 
                fopen(port);
                if isempty(CommPort)
                    CommPort = port;
                    CommPortIndex(i) = 1;
                else CommPort(end+1) = port;
                    CommPortIndex(i) = length(CommPort);
                end
            end
            for j=1:1:4
                a = Plant.Generator(i).VariableStruct.Measure.(char(Names(j)));
                if a~=0 
                    port = udp('192.168.0.8', a, 'LocalPort',a+10000); %don't open
                    if isempty(MeasurePort)
                        MeasurePort = port;
                        MeasurePortIndex(i,j) = 1;
                    else MeasurePort(end+1) = port;
                        MeasurePortIndex(i,j) = length(MeasurePort);
                    end
                end
            end
        end
    end     
end