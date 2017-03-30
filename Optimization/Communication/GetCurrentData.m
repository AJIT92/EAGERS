function RealData = GetCurrentData(Date)
global RealTimeData RealTime Virtual DispatchWaitbar
if Date ~=inf && Date>(RealTimeData.Timestamp(end-1))% End Operation or Simulation
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
    x1 = max(1,nnz(RealTimeData.Timestamp<=Date));
    S = fieldnames(RealTimeData.Demand);
    RealData.Timestamp = RealTimeData.Timestamp(x1);
    RealData.Temperature = RealTimeData.Temperature(x1);
    for i = 1:1:length(S)
        RealData.Demand.(S{i}) = RealTimeData.Demand.(S{i})(x1,:);
    end
end