classdef plant
    %both planning and dispatch tools use plant objects named Plant
    properties
        Name = 'New Plant';
        Data = []; %stores historical profiles and planning tool results 
        Generator = [];
        Dispatch = []; %recording of optimized dispatch
        optimoptions = []; %options for optimization times/CHP/demands etc.
        OpMatA = []; %optimization matrices for fit A
        OpMatB = []; %optimization matrices for fit B
        Threshold = []; %thresholds for when a generator needs to turn off
        OneStep = []; %optimization for one step (used for unit commitment solving)
        Online = []; %optimization for online control
        Control = [];
        Economic = [];
        Plotting = [];%these are color maps
        GUIhandles = [];
        Network = [];
        nodalVar = [];
    end
    methods 
        function obj = plant(name)
            obj.Name = name;
            %add data structure
            obj.Data.Timestamp = [];
            obj.Data.Holidays = [];
            obj.Data.Temperature = [];
            obj.Data.NameList = [];
            obj.Data.Demand = [];
            obj.Data.Demand.E = [];
            obj.Data.PlanningResult = [];
            obj.Data.HistProf = [];
            obj.Data.Weather = [];
            obj.Data.PlanningResult = [];
            %add generators
            obj.Generator = [];
            obj.Generator.Type = [];
            obj.Generator.Name = [];
            obj.Generator.Source = [];
            obj.Generator.Output = [];
            obj.Generator.Size = [];
            obj.Generator.VariableStruct = [];
            obj.Generator.OpMatA = [];
            obj.Generator.OpMatB = [];
            obj.Generator.Enabled = [];
            %add dispatch structure
            obj.Dispatch = [];
            obj.Dispatch.Dispatch = [];
            obj.Dispatch.RunData = [];
            obj.Dispatch.Predicted = [];
            obj.Dispatch.Economic = [];
            obj.Dispatch.ActiveGenerators = {};
            %add optimization options
            obj.optimoptions.Interval = 1;
            obj.optimoptions.Horizon = 24;
            obj.optimoptions.Resolution = 1;
            obj.optimoptions.Topt = 3600;
            obj.optimoptions.Tmpc = 600;
            obj.optimoptions.nsSmooth = 0;
            obj.optimoptions.scaletime = 1;
            obj.optimoptions.fastsimulation = 1;
            obj.optimoptions.tspacing = 'constant';
            obj.optimoptions.sequential = 1;
            obj.optimoptions.excessHeat = 1;
            obj.optimoptions.Outputs = {'E'};
            obj.optimoptions.scaleTime = 1;
            obj.optimoptions.thresholdSteps = 6;
            obj.optimoptions.Buffer = 20;
            %add control options for planning tool
            obj.Control.Name = '3_DiurnalPeak';
            obj.Control.Type = 'Control';
            %add economic characteristics for planning tool
            obj.Economic.InstallCost = 0;
            obj.Economic.Incentive = 0;
            obj.Economic.StackReplaceCost = 0;
            obj.Economic.FinanceYrs = 1;
            obj.Economic.DepreciateYrs = 5;
            obj.Economic.Inflation = 1;
            obj.Economic.Interest = 3;
            obj.Economic.Payback = 5;
            obj.Economic.LifeYrs = 20;
            obj.Economic.LifekWh = 20000;
            obj.Economic.ProdCosts = [];            
        end
        function updatedPlant = updatePlantVersion(Plant)
        end
    end
end
        