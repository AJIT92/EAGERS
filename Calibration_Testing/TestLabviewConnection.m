port1 = udp('192.168.0.8', 31100, 'LocalPort',41100); %mGT power [kW]
port2 = udp('192.168.0.8', 31200, 'LocalPort',41200); %mGT fuel
port3 = udp('192.168.0.8', 31400, 'LocalPort',41400); %ICE power [kW]
port4 = udp('192.168.0.8', 32100, 'LocalPort',42100); %SOC [kJ]
port5 = udp('192.168.0.8', 31500, 'LocalPort',41500); %ICE fuel [kW]
port9 = udp('192.168.0.8', 31900, 'LocalPort',41900); %fan thermal power [kW]

for k = 1:1:10
    tic
    T = timer('TimerFcn',@timerFCN,'StartDelay',3);
    start(T)
    
    fopen(port1);
    fopen(port2);
    fopen(port3);
    fopen(port4);
    fopen(port5);
    fopen(port9);
    mGT =(str2double(char(fread(port1))));
    fuel1 =(str2double(char(fread(port2))));
    ICE =(str2double(char(fread(port3))));
    SOC =(str2double(char(fread(port4))));
    fuel2 =(str2double(char(fread(port5))));
    Fan =(str2double(char(fread(port9))))
    fclose(port1);
    fclose(port2);
    fclose(port3);
    fclose(port4);
    fclose(port5);
    fclose(port9);
    wait(T)
    delete(T)
    toc
end

delete(port9);
delete(port1);
delete(port3);
delete(port2);
delete(port5);
delete(port4);


clear port9;
clear port1;
clear port3;
clear port2;
clear port5;
clear port4;