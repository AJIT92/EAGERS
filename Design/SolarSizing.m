function SolarSizing(Index,ScaleTo,Repeat,NewSize)
global Project

AnnualPower = zeros(length(Project.Renewable.Solar),1);
for i = 1:1:length(Project.Renewable.Solar)
    solar = Project.Renewable.Solar(i);
    if strcmp(solar.Tracking,'fixed')
        AnnualPower(i) = 8760*solar.Sizem2*sum((solar.Irrad/1000).*cosd(solar.SunZen-solar.Tilt).*(cosd(solar.SunAz-solar.Azimuth))*solar.Eff)/length(solar.SunAz);
    elseif strcmp(solar.Tracking,'1axis')
        AnnualPower(i) = 8760*solar.Sizem2*sum((solar.Irrad/1000).*cosd(solar.SunZen-solar.Tilt)*solar.Eff)/length(solar.SunAz);
    else AnnualPower(i) = 8760*solar.Sizem2*sum((solar.Irrad/1000)*solar.Eff)/length(solar.SunAz);
    end
end
solar = Project.Renewable.Solar;
if strcmp(ScaleTo,'demand')
    Annual = 8760*sum(Project.Building.DemandE)/length(Project.Building.DemandE);
elseif strcmp(ScaleTo,'generation')
    CHPsize = 0;
    for i = 1:1:length(Project.System.CHP)
        CHPsize = CHPsize + Project.System.CHP(i).SysSize(1);
    end
    Annual = 8760*CHPsize;
end
if strcmp(Repeat,'all')
    oldSize = sum(AnnualPower)/Annual*100;
else oldSize = AnnualPower(Index)/Annual*100;
end
if oldSize ==0
    scale = 0;
else scale = (NewSize/oldSize);
end
if strcmp(Repeat,'all')
    for i = 1:1:length(solar)
        solar(i).Size = solar(i).Size*scale;
        solar(i).Sizem2= solar(i).Sizem2*scale;
        solar(i).DemFrac = solar(i).DemFrac*scale;
        solar(i).GenFrac = solar(i).GenFrac*scale;
    end
else
    solar(Index).Size = solar(Index).Size*scale;
        solar(Index).Sizem2= solar(Index).Sizem2*scale;
        solar(Index).DemFrac = solar(Index).DemFrac*scale;
        solar(Index).GenFrac = solar(Index).GenFrac*scale;
end
Project.Renewable.Solar = solar;