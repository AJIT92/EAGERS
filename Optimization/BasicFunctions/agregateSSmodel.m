%% Build agregate state space model:
function agregateSSmodel(SSi)
%SSmpc: the agregated state-space model seen by the MPC, the time may be scaled by both Tmpc and scaletime.
%GenSSindex: A matrix used to relate the primary and secondary states of each generator to their index in the states of the agregated state space model SSmpc.
%GendX: The derivative of the generator state, assume steady initial condition and GendX=0. Need to keep track of this for MPC
global Plant GenSSindex SSmpc GendX
nG = length(Plant.Generator);
GenSSindex = [];
gen = [];
nCHP = 0;
for i=1:1:nG
    if ~isempty(SSi) && ~isempty(SSi(i).A)
        gen(end+1) = i;
        if isfield(Plant.Generator(i).OpMatB.output,'H') && isfield(Plant.Generator(i).OpMatB.output,'E')
            nCHP = nCHP+1;
        end
    end
end
GenSSindex.Primary = zeros(length(gen),2);
GenSSindex.Secondary = zeros(nCHP,2);
A = [];
B = [];
C = [];
r = 0;
w = 0;
chp=0;
for i = 1:1:length(gen)
    sA = size(SSi(gen(i)).A);
    if sA(1) ==2 %SISO 2nd order response (1 output)
        A = [A zeros(w,2); zeros(2,w) SSi(gen(i)).A];
        B = [B zeros(w,1); zeros(2,r) SSi(gen(i)).B];
        C = [C zeros(w/2,2); zeros(1,w) SSi(gen(i)).C];
        GenSSindex.Primary(i,1) = w+1; %1st column gives index in state vector x
        GenSSindex.Primary(i,2) = w/2+1; %2nd column gives index in output (y = Cx)
        w = w+2; %total states of X
    elseif sA(1) ==4 %SIMO 2nd order response (2 outputs)
        chp = chp+1;
        A = [A zeros(w,4); zeros(4,w) SSi(gen(i)).A];
        B = [B zeros(w,1); zeros(4,r) SSi(gen(i)).B];
        C = [C zeros(w/2,4); zeros(2,w) SSi(gen(i)).C];
        GenSSindex.Primary(i,1) = w+1; %1st column gives index in state vector x
        GenSSindex.Primary(i,2) = w/2+1; %2nd column gives index in output (y = Cx)
        GenSSindex.Secondary(chp,1) = w+3;  %1st column gives index in state vector x
        GenSSindex.Secondary(chp,2) = w/2+2; %2nd column gives index in output (y = Cx)
        w = w+4;%total states of X
    end 
    r = r+1; %input #
end
if ~isempty(SSi)
    SS.A = A;
    SS.B = B;
    SS.C = C;
    SS.D = zeros(length(C(:,1)),length(B(1,:)));
    SS.Dt = 1; %sampling time
    % SSmpc = changeTimestep(SS,Plant.optimoptions.Tmpc,1);
    % GendX = zeros(1,nG);
end