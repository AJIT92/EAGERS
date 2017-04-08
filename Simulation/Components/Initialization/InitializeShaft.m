function block = InitializeShaft(varargin)
% a simple shaft model
% Three (3) inlets: WTurbine, WCompressor, and Gen_Power
% One (1) outlets: RPM, Steady_Power
% One (1) states: Shaft Speed
global modelParam;
block = varargin{1};
if length(varargin)==1 % first initialization

    block.PMoI = pi()/2*block.Radius^4;%Shaft Polar Moment of Inertia, relates shaft speed to work
    
    BlockPort1 = block.connections{1};
    if ~isempty(BlockPort1) && ~strcmp('Tags',BlockPort1(1:4))
        r = strfind(BlockPort1,'.');
        connectedBlock1 = BlockPort1(1:r-1);
    end
    if isfield(modelParam,connectedBlock1) && isfield(modelParam.(connectedBlock1),'RPMdesign')
        block.Scale = modelParam.(connectedBlock1).RPMdesign/60*(2*pi); %shaft speed in Rad/s normalized by design shaft speed.
    else
        BlockPort2 = block.connections{2};
        if ~isempty(BlockPort2) && ~strcmp('Tags',BlockPort2(1:4))
            r = strfind(BlockPort2,'.');
            connectedBlock2 = BlockPort2(1:r-1);
        end
        if  isfield(modelParam,connectedBlock2) && isfield(modelParam.(connectedBlock2),'RPMdesign')
            block.Scale = modelParam.(connectedBlock2).RPMdesign/60*(2*pi); %shaft speed in Rad/s
        else block.Scale = (block.RPMinit/60*(2*pi));
        end
    end
    block.IC = (block.RPMinit/60*(2*pi))/block.Scale; %shaft speed in Rad/s
    
    block.InletPorts = {'WTurbine','WCompressor','Gen_Power'};
    block.WTurbine.IC = 200;%in KW
    block.WCompressor.IC = 100;%in KW
    block.Gen_Power.IC = 100;%in KW

    block.OutletPorts = {'RPM','Steady_Power'};
    block.RPM.IC = block.RPMinit; %shaft speed in RPM
    block.Steady_Power.IC = 100;%in KW
    
    block.P_Difference = {};
    %no dMdP or mFlow (no pressure calculations)
end
if length(varargin)==2 %% Have inlets connected, re-initialize
    Inlet = varargin{2};
    block.Steady_Power.IC = Inlet.WTurbine - Inlet.WCompressor;%in KW
end