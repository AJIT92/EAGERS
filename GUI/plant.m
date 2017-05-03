classdef plant
    %both planning and dispatch tools use plant objects named Plant
    properties
        Name = 'New Plant';
        Data = []; %stores historical profiles and planning tool results 
        optimoptions = []; %options for optimization times/CHP/demands etc.
        Generator = [];
        Network = [];
%         Plotting = [];%these are color maps
%         GUIhandles = [];
%         Dispatch = []; %recording of optimized dispatch
%         OpMatA = []; %optimization matrices for fit A
%         OpMatB = []; %optimization matrices for fit B
%         Threshold = []; %thresholds for when a generator needs to turn off
%         OneStep = []; %optimization for one step (used for unit commitment solving)
%         Online = []; %optimization for online control
    end
    methods 
        function obj = plant(name)
            obj.Name = name;
            %add data structure
            obj.Data.Timestamp = [];
            obj.Data.Holidays = [];
            obj.Data.Temperature = [];
            obj.Data.Demand = [];
            obj.Data.HistProf = [];
%             obj.Data.PlanningResult = [];
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
            obj.optimoptions.sequential = 0;
            obj.optimoptions.excessHeat = 1;
            obj.optimoptions.scaleTime = 1;
            obj.optimoptions.thresholdSteps = 6;
            obj.optimoptions.Buffer = 20;
            obj.optimoptions.method = 'Dispatch';
            obj.optimoptions.MixedInteger = true;
            obj.optimoptions.SpinReserve = false;
            obj.optimoptions.SpinReservePerc = 0;
            obj.optimoptions.Outputs = {};
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
            %add Network structure
            obj.Network.name = '';
            obj.Network.Equipment = {};
            
        end
        function TrainHistoricalData
            %update historical surface fits
        end
    end
end
        