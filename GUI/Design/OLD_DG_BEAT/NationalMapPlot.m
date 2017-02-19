function NationalMapPlot(Y,states2,units,ColorScaleMax,ColorScaleMin,ticks,color)
global Model_dir MapHandle
figure(MapHandle)
clf(MapHandle)
mapAx = USA48map(states2);
Cmap = colormap(color);

set(mapAx, 'Visible', 'off','Fontsize',14)
gridm('off')
latlim = getm(mapAx, 'MapLatLimit');
lonlim = getm(mapAx, 'MapLonLimit');
dataFileName = fullfile(Model_dir, 'GUI','Mapping','MapData','usastatehi'); 
states = shaperead(dataFileName,...
    'UseGeoCoords', true, 'BoundingBox', [lonlim', latlim']);

states3 = [];
for j = 1:1:length(states)
    if max(strcmp(states(j).Name,states2))==1
        states3(end+1) = j;
    end
end
for i = 1:1:length(Y)
    StateColor = Cmap(max(1,floor(length(Cmap)*((Y(i)-ColorScaleMin)/(ColorScaleMax-ColorScaleMin)))),1:3);
    geoshow(mapAx, states(states3(i)), 'FaceColor', StateColor)
end
% h = colorbar('CLim',[ColorScaleMin ColorScaleMax],'YTick',linspace(0,1,ticks),'YTickLabel',linspace(ColorScaleMin,ColorScaleMax,ticks),'FontSize',14);
h = colorbar('YTick',linspace(0,1,ticks),'YTickLabel',round(linspace(ColorScaleMin,ColorScaleMax,ticks)),'FontSize',14);
% caxis([ColorScaleMin ColorScaleMax]); % 2014b recode
ylabel(h,units,'FontSize',16);
%% http://www.mathworks.com/help/matlab/ref/colorbar.html

function h = USA48map(states)
%REGIONMAP Construct a map axes for a region of the USA.
global Model_dir
mapDataDir = fullfile(Model_dir,'GUI','Mapping','MapData'); 
load(fullfile(mapDataDir,'regions'),'stateRegions');
validRegions =  stateRegions;
mapLimitIncrement = 1;  % Degrees

% Set the ellipsoid to a reference sphere that models the Earth.
e = referenceSphere('earth','meters');

% Lookup the latitude/longitude limits for the region and
% combine limits from multiple regions.
[latlim, lonlim] = getLimitsFromRegions(...
    states, validRegions, mapLimitIncrement);
h = constructMapAxesUSA(latlim, lonlim, e);     

%--------------------------------------------------------------------------

function [latlim, lonlim] = getLimitsFromRegions(...
    regionNames, validRegions, inc)
% Return the most compact possible latitude and longitude limits that
% encompass the regions listed by name in cell array regionNames.

latlim = [];
lonlim = [];

validRegionNames = {validRegions.name};

for k = 1:numel(regionNames)
    index = strncmpi(regionNames{k},validRegionNames,numel(regionNames{k}));
    if ~any(index)
        error(message('map:regionmap:unknownRegion',regionNames{k}))
    end
    index = find(index);
    index = index(1);   % Should be a scalar, but just in case...
    latlim = mergelatlimits(latlim, validRegions(index).latlim);
    lonlim = mergelonlimits(lonlim, validRegions(index).lonlim);
    % Possible enhancements:  Merge all the regions simultaneously to
    % ensure that the longitude limit result is independent of the
    % processing order.  Also ensure that the longitude limits do not span
    % more than 360 degrees.
end

% Snap map limits to increments of INC, with a 1 degree buffer, except
% for limits that are already exact multiples of INC.
buffer = 1;

if mod(latlim(1),inc) ~= 0
    latlim(1) = inc * floor((latlim(1) - buffer)/inc);
end
if mod(lonlim(1),inc) ~= 0
    lonlim(1) = inc * floor((lonlim(1) - buffer)/inc);
end
if mod(latlim(2),inc) ~= 0
    latlim(2) = inc * ceil((latlim(2) + buffer)/inc);
end
if mod(lonlim(2),inc) ~= 0
    lonlim(2) = inc * ceil((lonlim(2) + buffer)/inc);
end

if latlim(1) < -90
    latlim(1) = -90;
end
if latlim(2) > 90
    latlim(2) = 90;
end

%--------------------------------------------------------------------------

function latlim = mergelatlimits(latlim1, latlim2)

% Compute the tightest possible latitude limits encompassing both the
% interval defined by 1-by-2 vector LATLIM1 and the interval defined by
% 1-by-2 vector in LATLIM2.  Note that either input could be empty.
limits = [latlim1 latlim2];
latlim = [min(limits) max(limits)];

%--------------------------------------------------------------------------

function lonlim = mergelonlimits(lonlim1, lonlim2)

% Compute the tightest possible longitude limits encompassing both the
% interval defined by 1-by-2 vector LONLIM1 and the interval defined by
% 1-by-2 vector LONLIM2.  In addition, LONLIM1, LONLIM2, or both may be
% empty.

if isempty(lonlim1)
    lonlim = lonlim2;  
elseif isempty(lonlim2)
   lonlim = lonlim1;
else
    % Shift both intervals such that the first one starts at zero.  Call
    % the shifted versions i1 and i2.
    s1 = lonlim1(1);
    i1 = zero22pi(lonlim1 - s1);
    i2 = lonlim2 - s1;
    if zero22pi(i2(1)) <= i1(2)
        % We have overlap, with interval 2 starting within interval 1
        s2 = i2(1) - zero22pi(i2(1));
        % If necessary, shift i2 by a multiple of 360 degrees, so that
        % i2(1) falls within i1.  Call the result j2.
        j2 = i2 - s2;
        % Merge i1 and j2
        j = [0, max(i1(2),j2(2))];
    elseif zero22pi(i2(2)) <= i1(2)
        % We have overlap, with interval 2 ending within interval 1
        s2 = i2(2) - zero22pi(i2(2));
        % If necessary, shift i2 by a multiple of 360 degrees, so that
        % i2(2) falls within i1.  Call the result j2.
        j2 = i2 - s2;
        % Merge i1 and j2
        j = [min(0,j2(1)) i1(2)];
    else
        % Neither overlap condition was met; there is no overlap. We can
        % define j (shifted output interval) by either putting i2 to the
        % east of i1, or by putting it to the west.  We'll make the choice
        % that minimizes the width of j.
        width1 = zero22pi(i2(2));   % Width putting i2 to the east
        width2 = i1(2) - (zero22pi(i2(1)) - 360);  % Width putting i2 to the west
        if width1 <= width2
            j = [0 width1];
        else
            j = i1(2) + [-width2 0];
        end
    end
    % Undo the shift s1
    lonlim = j + s1;
end

% The following changes the behavior of usamap('alaska') and seems
% unnecessary:
%
% lonlim = npi2pi(lonlim);  Return results in the interval [-180 180].

%--------------------------------------------------------------------------

function setCommonMapAxesProperties(ax, ticksize, roundat)
% Set common map axes properties.

setm(ax, ...
    'Frame', 'on',...
    'Grid',  'on',...
    'LabelRotation', 'on',...
    'MeridianLabel', 'on',...
    'ParallelLabel', 'on',...
    'MLabelParallel', 0, ...
    'MLabelLoc', ticksize,...
    'MLineLoc',  ticksize,...
    'PLabelLoc', ticksize,...
    'PLineLoc',  ticksize,...
    'MLabelRound', roundat,...
    'PLabelRound', roundat,...
    'GColor', [.75 .75 .75],...
    'GLineStyle', ':');

%--------------------------------------------------------------------------

function h = constructMapAxesUSA(latlim, lonlim, e)
% Construct a map axes suitable to the specified latitude and longitude
% limits, for maps covering all or part of the Conterminous U.S.

% Ensure row vectors.
latlim = latlim(:)';
lonlim = lonlim(:)';

% Ensure ascending lonlim.
if lonlim(1) > lonlim(2)
    if lonlim(1) > 180
        lonlim(1) = lonlim(1) - 360;
    else
        lonlim(2) = lonlim(2) + 360;
    end
end

% Compute a nice increment for labels and grid.
% Pick a ticksize that gives at least 3 grid lines

mindiff = min(abs([diff(latlim) diff(lonlim)]));
ticcandidates = [.1 .5 1 2 2.5 5 10 15 20 30 45 ] ;
[~,indx] = min( abs( 4 - (mindiff ./ ticcandidates) ));

ticksize = ticcandidates(indx);
roundat = 0;
if mod(ticksize,1) ~= 0; roundat = -1; end

% States are small enough to display conformally with little error
projection = 'lambert';

% set up the map axes
h = axesm(...
    'MapProjection',projection,...
    'MapLatLimit',latlim,...
    'MapLonLimit',lonlim,...
    'MapParallels',[], ...
    'Spheroid', e);
    
setCommonMapAxesProperties(h, ticksize, roundat);
   
set(h, 'Visible', 'off')
set(get(h,'Title'), 'Visible', 'on')
