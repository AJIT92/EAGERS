function block = InitializeControlSOFCsystem(varargin)
% Controls for SOFC system with external reformer, control heat exchanger bypass, blower power, anode recirculation and fuel flow rate
% Four (4) inlets: T hot, T cold, T average, Voltage
% Five (5) outlets:  Heater bypass, blower power, fuel flow rate, current, anode recirculation
% Three (3) states: Heater bypass, blower power, Current 
global F
block = varargin{1};
if length(varargin)==1 %first initialization
    block.IC = [1 1 1]; % inital condition
    block.Scale = block.IC;
    block.description = {'Heater Bypass';'Blower Power';'Fuel Cell Current';};
    
    Target = zeros(length(block.Target),1);
    for j = 1:1:length(Target)
        Target(j) = lookupVal(block.Target{j});
    end
    block.Target  = Target;
    
    block.Cells = lookupVal(block.Cells);
    block.Fuel = lookupVal(block.Fuel);
    block.Utilization = lookupVal(block.Utilization);
    Recirculation = lookupVal(block.InitialAnodeRecirc);
    
    for j = 1:1:length(block.Scale)
        block.Scale(j) = lookupVal(block.InitConditions{j});
    end
    initCurrent = block.Scale(3);
    
    block.IC = [block.Scale(1);1;1;]; % inital condition
    block.Scale = [1;block.Scale(2);block.Scale(3)];
    FuelFlow = (block.Cells*initCurrent/(2*F*block.Utilization*(4*block.Fuel.CH4+block.Fuel.CO+block.Fuel.H2))/1000);
    if Recirculation > 0
        %use reciculation from 1st initialization to set 
        Steam2Carbon = block.Target(3);
        if Steam2Carbon ==2 %steam to carbon ratio is unaffected by WGS effectiveness
            block.WGSeffective = 0.7;
        else
            eWGS = .7; %Initial guess of WGS efftiveness
            CH4 = block.Fuel.CH4*FuelFlow;
            r = Recirculation; 
            dWGS = 1e-6;
            error = 1;
            % Inlet = (Inlet + generated - consumed)*r  + New, thus inlet = New/(1-r) + (generated - consumed)*r/(1-r)
            while abs(error)>1e-6
                COin = block.Fuel.CO*FuelFlow/(1-r) + (block.Fuel.CH4 - eWGS*(block.Fuel.CH4+block.Fuel.CO))*FuelFlow*r/(1-r);
                Inlet.FuelMix.H2O = block.Fuel.H2O*FuelFlow/(1-r) + (block.Cells*initCurrent/(2*F*1000) - (block.Fuel.CH4 + (block.Fuel.CH4 + block.Fuel.CO)*eWGS)*FuelFlow)*r/(1-r);
                S2C = Inlet.FuelMix.H2O/(CH4 + 0.5*COin);
                error = Steam2Carbon - S2C;
                eWGS2 = eWGS+dWGS;
                COin2 = block.Fuel.CO*FuelFlow/(1-r) + (block.Fuel.CH4 - eWGS2*(block.Fuel.CH4+block.Fuel.CO))*FuelFlow*r/(1-r);
                H2Oin2 = block.Fuel.H2O*FuelFlow/(1-r) + (block.Cells*initCurrent/(2*F*1000) - (block.Fuel.CH4 + (block.Fuel.CH4 + block.Fuel.CO)*eWGS2)*FuelFlow)*r/(1-r);
                S2C2 = H2Oin2/(CH4 + 0.5*COin2);
                dSdWGS = (S2C2 - S2C)/dWGS;
                if dSdWGS ==0
                    error = 0;
                    eWGS = 1;
                else
                    eWGS = max(0,eWGS + error/dSdWGS);
                end
            end
            block.WGSeffective = eWGS;
        end
    end
    
    block.PortNames = {'Hot','Cold','PEN_Temp','Voltage','HeaterBypass','Blower','Current','AnodeRecirc','FuelFlow'};
    block.Hot.type = 'in';
    block.Hot.IC = block.Target(2)+.5*block.Target(1); 
    block.Cold.type = 'in';
    block.Cold.IC = block.Target(2)-.5*block.Target(1); 
    block.PEN_Temp.type = 'in';
    block.Voltage.type = 'in';
    block.Voltage.IC = 0.85; 
    block.HeaterBypass.type = 'out';
    block.HeaterBypass.IC = block.IC(1)*block.Scale(1);
    block.Blower.type = 'out';
    block.Blower.IC = block.IC(2)*block.Scale(2); 
    block.Current.type = 'out';
    block.Current.IC =  block.IC(3)*block.Scale(3);
    block.AnodeRecirc.type = 'out';
    block.AnodeRecirc.IC = Recirculation;
    block.FuelFlow.type = 'out';
    block.FuelFlow.IC = FuelFlow;
    
    block.P_Difference = {};
    
    for i = 1:1:length(block.PortNames)
        if length(block.connections)<i || isempty(block.connections{i})
            block.(block.PortNames{i}).connected={};
        else
            if ischar(block.connections{i})
                block.(block.PortNames{i}).connected = block.connections(i);
            else
                block.(block.PortNames{i}).IC = block.connections{i};
                block.(block.PortNames{i}).connected={};
            end
        end
    end
%     blowPow = block.Blower.IC
elseif length(varargin)==2 %% Have inlets connected, re-initialize
    Inlet = varargin{2};
    
    NetPower = PowerDemandLookup(0);
    PowerError = (NetPower-block.Current.IC*Inlet.Voltage*block.Cells/1000)/NetPower;
    block.Current.IC = block.Current.IC*(1 + PowerError);
    
    averageT = mean(Inlet.PEN_Temp); %need to adjust either target temp of heat exchanger or deltaT target if heat exchanger is sizing automatically to hit a target temp
    TavgError = (block.Target(1) - averageT)/block.Target(2);
%     if Inlet.Cold<(block.Target(1) - 1.5*block.Target(2))
%         Inlet.Cold  = block.Target(1) -.75*block.Target(2);
%     end
    deltaT = (mean(Inlet.Hot)-mean(Inlet.Cold));
    dTerror =(deltaT-block.Target(2))/block.Target(2);

    Steam2Carbon = block.Target(3);
%     dTerror
%     TavgError
%     if abs(dTerror)>.1
%         block.Blower.IC = block.Blower.IC*(1+.5*dTerror); % change blower power to change flow rate
%     else block.Blower.IC = block.Blower.IC*(1+dTerror)^2; % change blower power to change flow rate
%     end
%     blowPow = block.Blower.IC
    
    %%if too cold, add more fuel, it will burn and more heat will be recovered
    Utilization = block.Utilization*(1-0.05*TavgError)
    block.Utilization = Utilization;
%     block.HeaterBypass.IC = block.HeaterBypass.IC;
    
    FuelFlow = (block.Cells*block.Current.IC/(2*F*block.Utilization*(4*block.Fuel.CH4+block.Fuel.CO+block.Fuel.H2))/1000);
    block.FuelFlow.IC = FuelFlow;
    
    if block.AnodeRecirc.IC > 0
        r = block.AnodeRecirc.IC; %Initial guess of anode recirculation
        CH4 = block.Fuel.CH4*FuelFlow;
        dr = 1e-6;
        error = 1;
        % Inlet = (Inlet + generated - consumed)*r  + New, thus inlet = New/(1-r) + (generated - consumed)*r/(1-r)
        while abs(error)>1e-6
            COin = block.Fuel.CO*FuelFlow/(1-r) + (block.Fuel.CH4 - block.WGSeffective*(block.Fuel.CH4+block.Fuel.CO))*FuelFlow*r/(1-r);
            Inlet.FuelMix.H2O = block.Fuel.H2O*FuelFlow/(1-r) + (block.Cells*block.Current.IC /(2*F*1000) - (block.Fuel.CH4 + (block.Fuel.CH4 + block.Fuel.CO)*block.WGSeffective)*FuelFlow)*r/(1-r);
            S2C = Inlet.FuelMix.H2O/(CH4 + 0.5*COin);
            error = Steam2Carbon - S2C;
            r2 = r+dr;
            COin2 = block.Fuel.CO*FuelFlow/(1-r2) + (block.Fuel.CH4 - block.WGSeffective*(block.Fuel.CH4+block.Fuel.CO))*FuelFlow*r2/(1-r2);
            H2Oin2 = block.Fuel.H2O*FuelFlow/(1-r2) + (block.Cells*block.Current.IC /(2*F*1000) - (block.Fuel.CH4 + (block.Fuel.CH4 + block.Fuel.CO)*block.WGSeffective)*FuelFlow)*r2/(1-r2);
            S2C2 = H2Oin2/(CH4 + 0.5*COin2);
            dSdr = (S2C2 - S2C)/dr;
            r = r + error/dSdr;
        end
        block.AnodeRecirc.IC = r;
    end
%     block.InitializeError = max(abs(TavgError),abs(dTerror)); 
    block.InitializeError = abs(TavgError);
    
    block.Scale = [1, block.Blower.IC, block.Current.IC];
    block.IC = [block.HeaterBypass.IC-TavgError*block.PropGain(1), 1-dTerror*block.PropGain(2),1-PowerError*block.PropGain(3)];
end

function const = lookupVal(initCond)
global modelParam
if ischar(initCond)
    r = strfind(initCond,'.');
    if ~isempty(r)
        r = [r,length(initCond)+1];
        A = modelParam.(initCond(1:r(1)-1));
        for i = 2:1:length(r)
            field = initCond(r(i-1)+1:r(i)-1);
            A = A.(field);
        end
    else
        A = modelParam.(initCond);
    end
    if isstruct(A) && isfield(A,'T')
        const = NetFlow(A);
    elseif length(A)>1
        const  = sum(A);
    else
        const = A;
    end
else const  = initCond;
end