function [GeneratorProfile] = PredictRenewable(D1,D2,Gen)
GeneratorProfile = zeros((D2-D1)*96+1,1);
A = datevec(D1);
B = datenum(A(1),1,1);
Xts = (D1-B)*96;
XFts = Xts+length(GeneratorProfile)-1;
%% Predict Renewable Generation
if Gen.Enabled == 1
    if strcmp(Gen.VariableStruct.Type,'Solar') 
        if strcmp(Gen.VariableStruct.Tracking,'fixed')
            Power = Gen.VariableStruct.Sizem2*(Gen.VariableStruct.Irrad/1000).*cosd(Gen.VariableStruct.SunZen-Gen.VariableStruct.Tilt).*(cosd(Gen.VariableStruct.SunAz-Gen.VariableStruct.Azimuth))*Gen.VariableStruct.Eff;
        elseif strcmp(Gen.VariableStruct.Tracking,'1axis')
            Power = Gen.VariableStruct.Sizem2*(Gen.VariableStruct.Irrad/1000).*cosd(Gen.VariableStruct.SunZen-Gen.VariableStruct.Tilt)*Gen.VariableStruct.Eff;
        else Power = Gen.VariableStruct.Sizem2*(Gen.VariableStruct.Irrad/1000)*Gen.VariableStruct.Eff; %Dual axis
        end

        GeneratorProfile = Power(Xts:XFts);
    elseif strcmp(Gen.VariableStruct.Type,'Wind')
%             GeneratorProfile = Power(Xts:XFts);
    end
end