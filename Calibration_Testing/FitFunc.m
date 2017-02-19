function [Data] = FitFunc( ~ )
global Plant Tags noModel

%%adjusts efficiencies and massflows to fit data

Options.ODESet = odeset('AbsTol',1e-6,'RelTol',1e-5);
%Options.ODESet = odeset('AbsTol',1e-5,'RelTol',1e-4);

Tags = [];
Tags.Options = [];
Tags.Options.ODESet = Options.ODESet;


SetPoint.MFlow = 0.49;
SetPoint.CompExitT = 482;
SetPoint.TET = 907.8;
SetPoint.PR = 4.26;
SetPoint.Exhaust = 582;
ErrorT = 0.001;
check = 0;
n = 1;
Val.temp = [1 1];

Val.MFlow = 0.4894;
Val.TurbEff = 0.7704;
Val.CompEff = 0.7626;
Val.HXEffectiveness = 0.8004;
Val.PR = 4.1212;
Val.IntGain = [0*4e-4; 1e-2; 4e-2;];
Val.PropGain = [0*8e-3; 5e-0; .75;];
noModel = 1;

MFlow(1) = Val.MFlow;
TurbEff(1) = Val.TurbEff;
CompEff(1) = Val.CompEff;
Plates(1) = Val.HXEffectiveness;
PR(1) = Val.PR;
while check == 0
n = n + 1;
    CompEff(n) = Val.CompEff;
    MFlow(n) = Val.MFlow;
    TurbEff(n) = Val.TurbEff;
    Plates(n) = Val.HXEffectiveness;
    PR(n) = Val.PR;

    run('GasTurbine2');
    
    if abs(Tags.Turb.MassFlow -SetPoint.MFlow)/SetPoint.MFlow >=ErrorT %Correct Massflow error
        MFlow(n) = Val.MFlow - 0.25*(Tags.Turb.MassFlow -SetPoint.MFlow)/SetPoint.MFlow;
    elseif abs(Tags.Comp.Temperature - SetPoint.CompExitT)/SetPoint.CompExitT >=ErrorT
        CompEff(n) = Val.CompEff + 0.9*(Tags.Comp.Temperature - SetPoint.CompExitT)/SetPoint.CompExitT;
    elseif abs(Tags.Turb.TET - SetPoint.TET)/SetPoint.TET >=ErrorT
        TurbEff(n) = Val.TurbEff + 0.25*(Tags.Turb.TET - SetPoint.TET)/SetPoint.TET;
    elseif abs(Tags.HX1.HotOut - SetPoint.Exhaust)/SetPoint.Exhaust >=ErrorT
        Plates(n) = Val.HXEffectiveness + 0.5*(Tags.HX1.HotOut - SetPoint.Exhaust)/SetPoint.Exhaust;
    elseif abs(Tags.Comp.PR - SetPoint.PR)/SetPoint.PR >=ErrorT
        PR(n) = Val.PR - 4.5*(Tags.Comp.PR - SetPoint.PR)/SetPoint.PR;
    else
        check = 1;
    end
    
    Val.MFlow = MFlow(n);
    Val.TurbEff = TurbEff(n);
    Val.CompEff = CompEff(n);
    Val.HXEffectiveness = Plates(n);
    Val.PR = PR(n);

    disp(Val.MFlow);
    disp(Val.TurbEff);
    disp(Val.CompEff);
    disp(Val.HXEffectiveness);
    disp(Val.PR);
end

noModel = 0;
Val.IntGain = [4e-4; 1e-2; 4e-2;];
Val.PropGain = [8e-3; 5e-0; .75;];
Val.temp = [1 0.2];
run('GasTurbine2');
global TagInf

Data(:,1) = TagInf.GTcontrol.GenPower;
Data(:,2) = TagInf.Shaft.RPM;
Data(:,3) = TagInf.Turb.MassFlow;
Data(:,4) = TagInf.GTcontrol.FuelFlow;
Data(:,5) = TagInf.GTcontrol.Efficiency;
Data(:,6) = TagInf.Comp.PR*101.25;
Data(:,7) = TagInf.Turb.TET;
Data(:,8) = TagInf.Comp.Beta;
Data(:,9) = TagInf.Comp.NRPM;
Data(:,10) = TagInf.Turb.Beta;
Data(:,11) = TagInf.Turb.NRPM;
Data(:,12) = Data(:,5)*100;
end