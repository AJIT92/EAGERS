function [R,Out] = KineticReformation(method,Flow,Pressure,KineticCoeff,Current,block)%find new reforming reaction rates
global Ru F
spec = fieldnames(Flow.Inlet);
spec = spec(~strcmp('T',spec));
K_WGS = exp(4189.8./Flow.Outlet.T -3.8242);% Water gas shift equilibrium constant
nout = NetFlow(Flow.Outlet);
% nin = NetFlow(Flow.Inlet);
% X_CH4 = (Flow.Inlet.CH4+Flow.Outlet.CH4)./(nin+nout)*Pressure*1000; %partial pressures in Pa
% X_H2O = (Flow.Inlet.H2O+Flow.Outlet.H2O)./(nin+nout)*Pressure*1000; %partial pressures in Pa
    X_CH4 = Flow.Outlet.CH4./nout*Pressure*1000; %partial pressures in Pa
    X_H2O = Flow.Outlet.H2O./nout*Pressure*1000; %partial pressures in Pa
if strcmp(method,'Achenbach')
    R.CH4 = KineticCoeff*(X_CH4.*exp(-8.2e4./(Ru*Flow.Outlet.T)));
elseif strcmp(method,'Leinfelder')
    R.CH4 = KineticCoeff*(30.8e10*X_CH4.*X_H2O.*exp(-2.05e5./(Ru*Flow.Outlet.T)));
elseif strcmp(method,'Drescher')
    R.CH4 = KineticCoeff*(288.52*X_CH4.*X_H2O.*exp(-1.1e4./(Ru*Flow.Outlet.T))./(1+16*X_CH4+0.143*X_H2O.*exp(3.9e4./(Ru*Flow.Outlet.T))));
elseif strcmp(method,'equilibrium')
    a = 4352.2./Flow.Outlet.T - 3.99;
    K_WGS = block.scaleK_WGS.*exp(a);% Water gas shift equilibrium constant
    K_CH4 = block.scaleK_CH4.*2459000.*exp(-6.187*a);
    CH4_eq = (Flow.Outlet.H2.^3.*Flow.Outlet.CO)./(K_CH4.*Flow.Outlet.H2O).*(Pressure./NetFlow(Flow.Outlet)).^2;
    R.CH4 = Flow.Inlet.CH4-CH4_eq;
end
%confirm we don't violate anything casuing negative species
R.CH4 = min([R.CH4,Flow.Inlet.CH4,Flow.Inlet.H2O],[],2);
R.CH4 = max([R.CH4,-Flow.Inlet.CO,-Flow.Inlet.H2/3],[],2);

CO_eq = Flow.Outlet.CO2.*Flow.Outlet.H2./(K_WGS.*Flow.Outlet.H2O);
R.WGS = (Flow.Inlet.CO+R.CH4)-CO_eq; %inlet CO + CO from reforming - outlet CO
R.WGS = min([R.WGS,Flow.Inlet.CO+R.CH4,Flow.Inlet.H2O-R.CH4],[],2);
R.WGS = max([R.WGS,-Flow.Inlet.CO2,-(Flow.Inlet.H2+3*R.CH4)],[],2);
%identify the target outflow that dY will converge to
for i = 1:1:length(spec)
    if strcmp(spec{i},'CO2') 
        Out.CO2 = Flow.Inlet.CO2 + R.WGS + Current.CO/(2*F*1000); 
        if strcmp(block.FCtype,'MCFC')
            Out.CO2 = Out.CO2 + (Current.H2+Current.CO)/(2*F*1000);
        end
    elseif strcmp(spec{i},'H2') 
        Out.H2 = (Flow.Inlet.H2 + 3*R.CH4 + R.WGS - Current.H2/(2*F*1000));
    elseif strcmp(spec{i},'H2O') 
        Out.H2O = (Flow.Inlet.H2O - R.CH4 - R.WGS + Current.H2/(2*F*1000)); 
    elseif strcmp(spec{i},'CH4') 
        Out.CH4 = (Flow.Inlet.CH4 - R.CH4); 
    elseif strcmp(spec{i},'CO') 
        Out.CO = (Flow.Inlet.CO + R.CH4 - R.WGS - Current.CO/(2*F*1000));  
    else
        Out.(spec{i}) = Flow.Inlet.(spec{i});  
    end
end 