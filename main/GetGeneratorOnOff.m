function [State] = GetGeneratorOnOff(Measure)
State = zeros(1,length(Measure));
for i = 1:1:length(State)
    if Measure(i) ~=0
        port = udp('192.168.0.8', Measure(i), 'LocalPort',Measure(i)+10000); %Pelgt generata [kW]
        fopen(port);
        State(i)=(str2double(char(fread(port))));
        fclose(port);
        delete(port);
    end
end