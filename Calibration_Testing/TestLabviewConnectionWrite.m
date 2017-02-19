port11 = udp('192.168.0.8', 51100, 'LocalPort',61100); %set MGT power (kW)
port12 = udp('192.168.0.8', 51200, 'LocalPort',61200); %set boiler thermal demand (kW)
port13 = udp('192.168.0.8', 51300, 'LocalPort',61300); %set ICE electrical demand (kW)
port15 = udp('192.168.0.8', 51500, 'LocalPort',61500); %set fan thermal demand (kW)
port17 = udp('192.168.0.8', 51700, 'LocalPort',61700); %mGT on
port18 = udp('192.168.0.8', 51800, 'LocalPort',61800); %ICE on

fopen(port11);
fopen(port12);
fopen(port13);
fopen(port15);
fopen(port17);
fopen(port18);


fwrite(port11,num2str(22),'char');%set MGT power (kW)
fwrite(port12,num2str(0),'char'); %set boiler thermal demand (kW)
fwrite(port13,num2str(7),'char'); %set ICE electrical demand (kW)
fwrite(port15,num2str(12),'char'); %set fan thermal demand (kW)

for k = 1:1:20
    tic
    T = timer('TimerFcn',@timerFCN,'StartDelay',1);
    start(T)
    wait(T)
    delete(T)

    fwrite(port11,num2str(20+k),'char');%set MGT power (kW)
    fwrite(port12,num2str(0+k),'char'); %set boiler thermal demand (kW)
    fwrite(port13,num2str(0+k),'char'); %set ICE electrical demand (kW)
    fwrite(port15,num2str(10+k),'char'); %set fan thermal demand (kW)
    k
    toc
end



fwrite(port17,num2str(1),'char'); %mGT on
fwrite(port18,num2str(0),'char'); %ICE on

fclose(port11);
delete(port11);
fclose(port12);
delete(port12);
fclose(port13);
delete(port13);
fclose(port15);
delete(port15);
fclose(port17);
delete(port17);
fclose(port18);
delete(port18);

clear port11;
clear port12;
clear port13;
clear port15;
clear port17;
clear port18;

% port4 = udp('192.168.0.8', 51400, 'LocalPort',61400); %??????
% fopen(port4);
% fwrite(port4,num2str(49),'char');

% port6 = udp('192.168.0.8', 51600, 'LocalPort',61600); %note electrical demand for log
% fopen(port6);
% fwrite(port6,num2str(72),'char');


% fclose(port4);
% delete(port4);

% fclose(port6);
% delete(port6);

