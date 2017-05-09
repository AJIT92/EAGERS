function RealData = GetCurrentData(Date)
global RealTimeData RealTime Virtual DispatchWaitbar
if Date ~=inf && any(Date>(RealTimeData.Timestamp(end-1)))% End Operation or Simulation
    if RealTime
        closePorts;
    end
    RealTime=0;%end condition for real simulation
    Virtual = 0;%end condition for virtual simulation
    close(DispatchWaitbar)
    T1 = timerfind('Name', 'dispTimer') ;
    T2 = timerfind('Name', 'optTimer') ;
    T3 = timerfind('Name', 'mpcTimer') ;
    T4 = timerfind('Name', 'fanTimer') ;
    Timers = [T1,T2,T3,T4];
    for i = 1:1:length(Timers)
        stop(Timers(i));
        delete(Timers(i))
    end
    RealData = [];
else 
    %% need better indexing system to call current data
    x1 = max(1,nnz(RealTimeData.Timestamp<Date(1))); 
    S = fieldnames(RealTimeData.Demand);
    RealData.Timestamp = Date;
    if length(Date) == 1
        r = (Date - RealTimeData.Timestamp(x1))/(RealTimeData.Timestamp(x1+1) - RealTimeData.Timestamp(x1));
        RealData.Temperature = (1-r)*RealTimeData.Temperature(x1) + r*RealTimeData.Temperature(x1+1);
        for i = 1:1:length(S)
            RealData.Demand.(S{i}) = (1-r)*RealTimeData.Demand.(S{i})(x1,:) + r*RealTimeData.Demand.(S{i})(x1+1,:);
        end
    else
        x2 = max(1,nnz(RealTimeData.Timestamp<Date(end))); 
        RealData.Temperature = interp1(RealTimeData.Timestamp(x1:x2+1),RealTimeData.Temperature(x1:x2+1),Date);
        for i = 1:1:length(S)
            RealData.Demand.(S{i}) = interp1(RealTimeData.Timestamp(x1:x2+1),RealTimeData.Demand.(S{i})(x1:x2+1,:),Date);
        end
    end
end