function SetGeneratorOnOff(newlyOn,newlyOff)
global CommPortOnOff CommPortOnOffIndex
for i = 1:1:length(newlyOn)
    if newlyOn(i) && CommPortOnOffIndex(i)~=0
        portNumber = CommPortOnOffIndex(i);
        port = CommPortOnOff(portNumber);
        fwrite(port,'1','char');
    end
end
for i = 1:1:length(newlyOff)
    if newlyOff(i) && CommPortOnOffIndex(i)~=0
        portNumber = CommPortOnOffIndex(i);
        port = CommPortOnOff(portNumber);
        fwrite(port,'0','char');
    end
end