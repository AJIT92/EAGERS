function Out = Shaft(t,Y, Inlet,block,string1)
% a simple shaft model
% Three (3) inlets: WTurbine, WCompressor, and Gen_Power
% Two (2) outlets: RPM, Steady_Power
% One (1) states: Shaft Speed
global Tags
if strcmp(string1,'Outlet')
    Out.Steady_Power = Inlet.WTurbine - Inlet.WCompressor;
    Out.RPM = Y(1)*60/(2*pi);%Converts from Radians per Second to RPM
    Tags.(block.name).RPM = Out.RPM;
elseif strcmp(string1,'dY')
    dY = 0*Y;
    ShaftPower = (Inlet.WTurbine - Inlet.WCompressor - Inlet.Gen_Power)*1000;%all units should be W
    Moment_Inertia = block.Density * block.Length * block.PMoI;%Moment of Inertia for the shaft
    dY(1) = ShaftPower/(Moment_Inertia*Y(1)); %dw/dt = P/(w*l) units of rad/s
    Out = dY;
end