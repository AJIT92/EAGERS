global Model_dir
Model_dir=strrep(which('EAGERS.m'),'EAGERS.m','');
p = path;
dir1 = fullfile(Model_dir,'Plant');
if isempty(strfind(p,dir1))% && isempty(strfind(p,'Plant'))
    SubPaths2Add = {fullfile(Model_dir);
                    fullfile(Model_dir,'data');
                    fullfile(Model_dir,'Calibration_Testing');
                    fullfile(Model_dir,'Design');
                    fullfile(Model_dir,'Design','ClimateZonesByState');
                    fullfile(Model_dir,'DesignProjects');
                    fullfile(Model_dir,'Emissions');
                    fullfile(Model_dir,'GUI');
                    fullfile(Model_dir,'GUI','Design');
                    fullfile(Model_dir,'GUI','Graphics');
                    fullfile(Model_dir,'GUI','HelpTool');
                    fullfile(Model_dir,'GUI','Mapping');
                    fullfile(Model_dir,'GUI','Optimization');
                    fullfile(Model_dir,'GUI','Simulation');
                    fullfile(Model_dir,'help');
                    fullfile(Model_dir,'Model Library');
                    fullfile(Model_dir,'Model Library','Saved Models');
                    fullfile(Model_dir,'Model Library','Saved Linearizations');
                    fullfile(Model_dir,'Optimization');
                    fullfile(Model_dir,'Optimization','BasicFunctions');
                    fullfile(Model_dir,'Optimization','ComplementaryQP');
                    fullfile(Model_dir,'Optimization','Communication');
                    fullfile(Model_dir,'Optimization','Forecasting');
                    fullfile(Model_dir,'Optimization','MPC');
                    fullfile(Model_dir,'Optimization','HydroQP');
                    fullfile(Model_dir,'Optimization','NetworkQP');
                    fullfile(Model_dir,'Optimization','NeuralNetwork');
                    fullfile(Model_dir,'Optimization','Threshold');
                    fullfile(Model_dir,'results');
                    fullfile(Model_dir,'results','LoggedData');
                    fullfile(Model_dir,'results','OptimalMap');
                    fullfile(Model_dir,'results','RealTime');
                    fullfile(Model_dir,'results','Virtualization');
                    fullfile(Model_dir,'Simulation');
                    fullfile(Model_dir,'Simulation','BasicFunctions');
                    fullfile(Model_dir,'Simulation','Components');
                    fullfile(Model_dir,'Simulation','Components','Initialization');
                    fullfile(Model_dir,'Simulation','CompressorMaps');
                    fullfile(Model_dir,'Simulation','Controls');
                    fullfile(Model_dir,'Simulation','Controls','Initialization');
                    fullfile(Model_dir,'Simulation','FCMaps');
                    fullfile(Model_dir,'Simulation','LookupFunctions');
                    fullfile(Model_dir,'System Library');
                    fullfile(Model_dir,'System Library','Boiler');
                    fullfile(Model_dir,'System Library','Buildings');
                    fullfile(Model_dir,'System Library','Buildings','RealBuildingData');
                    fullfile(Model_dir,'System Library','Chiller');
                    fullfile(Model_dir,'System Library','CHP Generator');
                    fullfile(Model_dir,'System Library','Control');
                    fullfile(Model_dir,'System Library','Economic');
                    fullfile(Model_dir,'System Library','Electric Generator');
                    fullfile(Model_dir,'System Library','Electric Storage');
                    fullfile(Model_dir,'System Library','Grid');
                    fullfile(Model_dir,'System Library','Heater');
                    fullfile(Model_dir,'System Library','HVAC');
                    fullfile(Model_dir,'System Library','NatGas');
                    fullfile(Model_dir,'System Library','Solar');
                    fullfile(Model_dir,'System Library','Thermal Storage');
                    fullfile(Model_dir,'System Library','Utility');
                    fullfile(Model_dir,'System Library','Weather'); 
                    fullfile(Model_dir,'System Library','Wind');
                    fullfile(Model_dir,'Plant','GenSets');
                    fullfile(Model_dir,'Plant','LoadProfiles');
                    fullfile(Model_dir,'Plant');};

    for i=1:length(SubPaths2Add)
        addpath(SubPaths2Add{i});
    end
end
close all
WelcomeScreen