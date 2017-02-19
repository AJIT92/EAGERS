function P = PowerDemandLookup(t)
global SimSettings
P = interp1(SimSettings.PowerTime,SimSettings.PowerDemand,t);
