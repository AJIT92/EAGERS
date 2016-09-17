function closePorts
global CommPort CommPortOnOff MeasurePort FanPortRead FanPortWrite AmbTempRead

if ~isempty(CommPortOnOff)
    for i =1:1:length(CommPortOnOff)
        port = CommPortOnOff(i);
        fclose(port);
        delete(port);
        clear port
    end
end
if ~isempty(CommPort)
    for i =1:1:length(CommPort)
        port = CommPort(i);
        fclose(port);
        delete(port);
        clear port
    end
end
if ~isempty(MeasurePort)
    for i =1:1:length(MeasurePort)
        port = MeasurePort(i);
%         fclose(port);
        delete(port);
        clear port
    end
end
if ~isempty(FanPortRead)
%     fclose(FanPortRead);
    delete(FanPortRead);
end
if ~isempty(AmbTempRead)
%     fclose(FanPortRead);
    delete(AmbTempRead);
end

if ~isempty(FanPortWrite)
    fclose(FanPortWrite);
    delete(FanPortWrite);
end
CommPort =[];
CommPortOnOff =[];
MeasurePort =[];
FanPortRead =[];
FanPortWrite =[];
AmbTempRead =[];