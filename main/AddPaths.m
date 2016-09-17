function AddPaths()
% This function adds paths needed to run EAGERS
global Model_dir 

%build up paths to add, start with subpaths using fullfile for proper
%directory syntax for multiple platforms
SubPaths2Add = {fullfile(Model_dir,'main');
                fullfile(Model_dir,'Plant');
                fullfile(Model_dir,'Plant','GenSets');
                fullfile(Model_dir,'Plant','LoadProfiles');
                fullfile(Model_dir,'component library');
                fullfile(Model_dir,'results');
                fullfile(Model_dir,'results','Virtualization');
                fullfile(Model_dir,'results','RealTime');
                fullfile(Model_dir,'help');
                fullfile(Model_dir,'data');
                fullfile(Model_dir,'component library','Boiler');
                fullfile(Model_dir,'component library','Buildings');
                fullfile(Model_dir,'component library','Chiller');
                fullfile(Model_dir,'component library','CHP Generator');
                fullfile(Model_dir,'component library','DispatchOptimizations');
                fullfile(Model_dir,'component library','Electric Generator');
                fullfile(Model_dir,'component library','Electric Storage');
                fullfile(Model_dir,'component library','Heater');
                fullfile(Model_dir,'component library','HVAC');
                fullfile(Model_dir,'component library','Solar');
                fullfile(Model_dir,'component library','Thermal Storage');
                fullfile(Model_dir,'component library','Utility');
                fullfile(Model_dir,'component library','Weather');
                fullfile(Model_dir,'component library','Wind');};

for i=1:length(SubPaths2Add)
    addpath(SubPaths2Add{i});
end

