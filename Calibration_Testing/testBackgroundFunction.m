function testBackgroundFunction

% Create a figure

figurehandle = figure                            ;
set(figurehandle, 'CloseRequestFcn', @myclosereq);
setappdata(0,'figurehandle',figurehandle)        ;

% Create initial data

delT = 1       ;%update frequency
X=0;
Y = 0;

setappdata(figurehandle,'delT', delT);

% Plot initial data to the figure

plot(X, Y);
xlabel('Time (s)')
ylabel('Cumulative of random generator')

% Create and start a timer

mytimer     = timer('TimerFcn'     ,@testBackground_update, ...
                   'StartDelay'   , 0            , ...
                   'Period'       , delT         , ...
                   'Name'         ,'objectTimer' , ...
                   'ExecutionMode','fixedrate'        );

start(mytimer);

% Disable access to the handle

set(figurehandle, 'HandleVisibility', 'callback')

function testBackground_update(varargin)
% Get handles

fig   = getappdata(0,'figurehandle');
axe   = get(fig,'Children');

% Find the axis data

h     = get(axe,'Children');

xdata = get(h,'XData');   
ydata = get(h,'YData');   

% Update figure for the new time

xdata = [xdata xdata(end)+1];
ydata = [ydata, ydata(end)+rand];

set(h, 'XData', xdata);
set(h, 'YData', ydata);

% Create a new figure title

title(axe,[ 'Cumulative of random'])

% Fix the 'tick spacing' problem when scrolling

set(axe, 'XLim', [min(xdata)  max(xdata)]);

% Make sure the axis limits go to zero

ylim    = get(axe,'YLim');
ylim(1) = 0              ;
ylim(2) = max(ydata)*1.2 ;
set(axe,'YLim',ylim)     ;
   

% --------------------------------------------------------------------------
% myclosereq()
% --------------------------------------------------------------------------
% Deletes the timer when the figure is closed 
% --------------------------------------------------------------------------
function myclosereq(varargin)

% Stop and delete timer

mytimer = timerfind('Name', 'objectTimer') ;

if ~isempty(mytimer)
  stop  (mytimer)                         ;
  delete(mytimer)                         ;
end

% Call normal close request

closereq                                   ;