function Flow = FCin2Out(T,Inlet,Dir, Type,Cells,Current,R,c_or_a) 
%flow of species out after steam reformation & water gas shift
%(does not normalize back to 1) this is important for finding G correctly)
global F
k = Dir(:,1);
r = length(k);
Flow.Outlet.T = T;
spec = fieldnames(Inlet);
if isempty(R)
    R.CH4 = zeros(length(T),1);
    R.WGS = zeros(length(T),1);
end
for i = 1:1:length(spec)
    if strcmp(spec{i},'T')
        Flow.Inlet.T(k,1) = Inlet.T;
    else
        Flow.Inlet.(spec{i})(k,1) = Inlet.(spec{i})/Cells/r;
    end
end
for j= 1:1:length(Dir(1,:))
    if j>1
        k2 = Dir(:,j);
        for i = 1:1:length(spec)
            if strcmp(spec{i},'T')
                Flow.Inlet.T(k2,1) = Flow.Outlet.T(k);
            else
                Flow.Inlet.(spec{i})(k2,1) = Flow.Outlet.(spec{i})(k);
            end
        end
        k = k2;
    end
    
    switch Type
        case {'SOFC';'MCFC'} %current is positive
            if strcmp(c_or_a,'anode') % H2 <-->H2O side (H2 consumption)
                side = 'H2';
            else side = 'O2';
            end
        case {'SOEC';'MCEC'} %current is negative
            if strcmp(c_or_a,'cathode') % H2 <-->H2O side (H2 production)
                side = 'H2';
            else side = 'O2';
            end
    end
    if strcmp(side,'H2')
        for i = 1:1:length(spec)
            if strcmp(spec{i},'CO2')
                switch Type
                    case 'SOFC'
                        Flow.Outlet.CO2(k,1) = Flow.Inlet.CO2(k)+R.WGS(k); %CO2 flow
                    case 'MCFC'
                        Flow.Outlet.CO2(k,1) = Flow.Inlet.CO2(k)+R.WGS(k) + Current(k)/(2*F*1000); %CO2 flow
                end
            elseif strcmp(spec{i},'CH4')
                Flow.Outlet.CH4(k,1) = Flow.Inlet.CH4(k)-R.CH4(k);%CH4 flow
            elseif strcmp(spec{i},'CO')
                Flow.Outlet.CO(k,1) = Flow.Inlet.CO(k)+R.CH4(k)-R.WGS(k); %CO flow
            elseif strcmp(spec{i},'H2')
                Flow.Outlet.H2(k,1) = Flow.Inlet.H2(k)+3*R.CH4(k)+R.WGS(k) - Current(k)/(2*F*1000); %H2 flow
            elseif strcmp(spec{i},'H2O')
                Flow.Outlet.H2O(k,1) = Flow.Inlet.H2O(k)-R.CH4(k)-R.WGS(k) + Current(k)/(2*F*1000); %H2O flow
            elseif ~strcmp(spec{i},'T')
                Flow.Outlet.(spec{i})(k,1) = Flow.Inlet.(spec{i})(k);
            end
            if abs(Flow.Outlet.(spec{i})(k,1))<1e-18 %zero
                Flow.Outlet.(spec{i})(k,1) = 0;
            end
        end
    elseif strcmp(side,'O2') 
        for i = 1:1:length(spec)
            if strcmp(spec{i},'CO2')
                switch Type
                    case 'SOFC'
                        Flow.Outlet.CO2(k,1) = Flow.Inlet.CO2(k); %CO2 flow
                    case 'MCFC'
                        Flow.Outlet.CO2(k,1) = Flow.Inlet.CO2(k) - Current(k)/(2*F*1000); %CO2 flow
                end
            elseif strcmp(spec{i},'O2')
                Flow.Outlet.O2(k,1) = Flow.Inlet.O2(k) - Current(k)/(4*F*1000); %O2 flow
            elseif ~strcmp(spec{i},'T')
                Flow.Outlet.(spec{i})(k,1) = Flow.Inlet.(spec{i})(k);
            end
            if abs(Flow.Outlet.(spec{i})(k,1))<1e-18 %zero
                Flow.Outlet.(spec{i})(k,1) = 0;
            end
        end
    end
end