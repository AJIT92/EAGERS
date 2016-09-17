function SetGeneratorOutput(GenSet)
global CommPort CommPortIndex
for i = 1:1:length(CommPortIndex)
    if CommPortIndex(i)~=0
        portNumber = CommPortIndex(i);
        port = CommPort(portNumber); 
        fwrite(port,num2str(GenSet(i)),'char');
    end
end