global Model_dir
Model_dir=strrep(which('CalibrateGT.m'),fullfile('BasicFunctions','CalibrateGT.m'),'');
addpath(Model_dir);
AddPaths()


RealData = xlsread('CapstoneTurnDownData');
RPMtol = 0.01;%error tolerance at RPM set points

check = 0;

%should calibrate RPM set points on comp/turb maps
while check == 0
    Data = FitFunc();
    turbine = load('RadialTurbineIterator');
    compressor = load('RadialCompressorIterator');
    TurbineRPM = turbine.map.RPM;
    CompressorRPM = compressor.map.RPM;

    turbine.map.RPM(turbine.map.RPM>1.00001) = [];
    turbine.map.RPM(turbine.map.RPM<min(Data(:,11))) = [];
    for i = 1:length(turbine.map.RPM)
        %find closest point to setpoint in turb map
        error = abs(Data(:,11) - turbine.map.RPM(i));
        [row,~] = find(error(:) == min(error(:)),1);
        %find point with same power output in RealData
        error = [];
        error = abs(RealData(:,1) - Data(row,1));
        [RealRow,~] = find(error(:) == min(error(:)),1);
        %find RPM error between model and real at set point
        RPMError(i) = (Data(row,2) - RealData(RealRow,2))/96000;%normalize to NRPM
        n = find(TurbineRPM == turbine.map.RPM(i),1);
        TurbineRPM(n) = TurbineRPM(n) - 0.2*RPMError(i);%move turb map RPM set point to minimize error
        if n > 1
            if TurbineRPM(n) < TurbineRPM(n-1)
                TurbineRPM(n-1) = TurbineRPM(n-1) - 0.2*RPMError(i);
            end
        end
    end
    RPMErrorT = abs(RPMError);
    error = [];
    RPMError = [];
    
    
    
    compressor.map.RPM(compressor.map.RPM>1.00001) = [];
    compressor.map.RPM(compressor.map.RPM<min(Data(:,9))) = [];
    for i = 1:length(compressor.map.RPM)
        %find closest point to setpoint in comp map
        error = abs(Data(:,9) - compressor.map.RPM(i));
        [row,~] = find(error(:) == min(error(:)),1);
        %find point with same power output in RealData
        error = [];
        error = abs(RealData(:,1) - Data(row,1));
        [RealRow,~] = find(error(:) == min(error(:)),1);
        %find RPM error between model and real at set point
        RPMError(i) = (Data(row,2) - RealData(RealRow,2))/96000;%normalize to NRPM
        n = find(CompressorRPM == compressor.map.RPM(i),1);
        CompressorRPM(n) = CompressorRPM(n) - 0.2*RPMError(i);%move comp map RPM set point to minimize error
        if n > 1
            if CompressorRPM(n) < CompressorRPM(n-1)
                CompressorRPM(n-1) = CompressorRPM(n-1) - 0.2*RPMError(i);
            end
        end
    end
    RPMErrorC = abs(RPMError);
    error = [];
    RPMError = [];
    
    
    if max(RPMErrorT) <= RPMtol && max(RPMErrorC) <= RPMtol
        check = 1;
    end
    turbine.map.RPM = TurbineRPM;
    compressor.map.RPM = CompressorRPM;
    map = turbine.map;
    pathname = strcat(Model_dir,'/CompressorMaps/');
    maploc = fullfile(pathname,'RadialTurbineIterator.mat');
    save(maploc,'map');
    map = [];
    map = compressor.map;
    maploc = fullfile(pathname,'RadialCompressorIterator.mat');
    save(maploc,'map');
    Data = [];
end