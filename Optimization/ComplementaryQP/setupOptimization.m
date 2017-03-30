%% build matrices for cQP
global Plant
Time = buildTimeVector(Plant.optimoptions);%% set up dt vector of time interval length
Plant.OpMatA = buildMatrices('A',Time); %build quadratic programming matrices for FitA
Plant.OpMatB = buildMatrices('B',Time);%build quadratic programming matrices for FitB
dt = Time - [0, Time(1:end-1)];
n = Plant.optimoptions.thresholdSteps;
Time2 = linspace(1,n,n)*dt(1)/n;
Plant.Threshold = buildMatrices('B',Time2);%build quadratic programming matrices used for threshold calculation with linear spacing from t =0 to t = dt(1).
Plant.OneStep = buildMatrices1Step(dt);%build quadratic programming matrices for 1Step at interval spacing of dt

A.Horizon = Time(1);%Plant.optimoptions.Resolution;%the horizon is not the resolution converted to seconds
A.Resolution = Plant.optimoptions.Topt/3600;%the resolution is the frequency of Topt
A.tspacing = 'constant';
Time3 = buildTimeVector(A);%% set up dt vector of time interval length
Plant.Online = [];
Plant.Online.QP = [];
Plant.Online.Organize = [];
Plant.Online.Timestamp = [];
for t = 1:1:length(Time3)
    Plant.Online(t) = buildMatrices('B',Time3(t:end)); %build the matrix for the onlineOptimLoop using FitB
end