function Out = ControlMCFCstack(t,Y, Inlet,block,string1)
% Controls for Fuel cell stack only, control air flow and inlet temperature and fuel flow rate
% Four (4) inlets: T hot, T cold, T average, Current
% Four (4) outlets: OxidentTemp, oxidant flow rate, fuel temp, fuel flow rate
% Two (2) states: OxidentTemp, oxidant flow rate, 
global Tags F
NetPower = PowerDemandLookup(t);
FuelFlow = block.Cells*sum(Inlet.Current)/(2*F*block.Utilization*(4*block.Fuel(1)+block.Fuel(2)+block.Fuel(4)))/1000;
averageT = mean(Inlet.PEN_Temp);
TavgError = (block.Target(1)-averageT)/block.Target(1);
deltaT = (mean(Inlet.Hot)-mean(Inlet.Cold));
dTerror =(deltaT-block.Target(2))/block.Target(2);
EquivError = (block.Target(3)-Inlet.Equiv);
if strcmp(string1,'Outlet')
    if t>0
        Out.OxidantTemp = Y(1)+(TavgError*block.PropGain(1))*block.Scale(1);
        Out.OxidantFlow = Y(2)+(dTerror*block.PropGain(2))*block.Scale(2);
        Out.FuelTemp = block.FuelTemp.IC;%currently no control of fuel inlet temp
        Out.FuelFlow = FuelFlow;
        Out.ValvePerc = Y(3)+(EquivError*block.PropGain(3))*block.Scale(3);
        Out.PowerStack = NetPower;
    else
        Out.OxidantTemp = block.OxidantTemp.IC;
        Out.OxidantFlow = block.OxidantFlow.IC; 
        Out.FuelTemp = block.FuelTemp.IC;
        Out.FuelFlow = block.FuelFlow.IC;
        Out.ValvePerc = block.ValvePerc.IC;
        Out.PowerStack = block.StackPower;
    end
elseif strcmp(string1,'dY')
    dY = Y*0;
    dY(1) = block.Gain(1)*TavgError;
    dY(2) = block.Gain(2)*dTerror;
    dY(3) = block.Gain(3)*EquivError;
    Out = dY;
    Tags.(block.name).Toxidant = Y(1)+(TavgError*block.PropGain(1))*block.Scale(1);
    Tags.(block.name).AirFlow = Y(2)+(dTerror*block.PropGain(2))*block.Scale(2);
    Tags.(block.name).Tfuel = block.FuelTemp.IC;%currently no control of fuel inlet temp
    Tags.(block.name).FuelFlow = FuelFlow;
    Tags.(block.name).ValvePerc = Y(3)+(EquivError*block.PropGain(3))*block.Scale(3);
end
