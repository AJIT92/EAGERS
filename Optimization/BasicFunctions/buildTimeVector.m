function Time = buildTimeVector(options)
%% set up dt vector of time interval length
if strcmp(options.tspacing,'constant') == 1
    nS = options.Horizon/options.Resolution; % # of time inervals
    Time = linspace(options.Resolution,options.Horizon,nS)';
%% if user manually specified timesteps
elseif strcmp(options.tspacing, 'manual')
    Time = options.manualT.*options.Horizon; %the user was told to specify steps as a portion of the horizon
    if options.manualT(end)~=1 %if they did not specify up to one horizon make the last entry one horizon
        Time(end+1) = options.Horizon;
    end
    roundTime = round(Time/(options.Tmpc/3600));%make sure that each step can contain a whole number of MPC steps
    Time = roundTime*(options.Tmpc/3600);
    [m,n] = size(Time);
    if m==1 & n>1
        Time = Time';
    end
%% if user specified logarithmic timesteps, then have more steps in the beggining, placing less certainty on the distant predictions
elseif strcmp(options.tspacing, 'logarithm')
    %the first step must be longer than the time for tmpc
    %made the first step a 5 min interval
    n1 = round(options.Horizon);
    if options.Topt/3600>1/12 %if the online loop is greater than 5 minutes, then make the shortest interval one onlineloop
        const = 1/(options.Topt/3600);
    else const = 12; %if the online loop runs more frequently than 5 min, then make the shortest interval 5 min
    end
    t = exp(1:n1)'/(const*exp(1));%divide by 12, because you want the first step to be 5 min
    Time = t(t<=24);
    roundTime = round(Time/(options.Tmpc/3600));%make sure that each step can contain a whole number of MPC steps
    Time = roundTime*(options.Tmpc/3600);
end
