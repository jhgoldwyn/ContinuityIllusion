function [t,x,I] = continuityModel(inputType,nTone,IT,IN,m,aE,aI,alph,bet,gon,goff,xinitial)

%  simulate x(t) firing rate dynamics for 
% "Dynamics of the Auditory Continuity Illusion"
% Cao, Parks, Goldwyn, 2021

% input parameters:
%   inputType: 'transient', 'sustained', or 'combined'.   or 'transient0' for special case of combined model parameter setting that receives only transient input
%   nTone: number of tones (0=noise only, 1=masking, 2=continuity)
%   IT: tone strength
%   IN: noise strength
%   aE: recurrent excitation
%   aI: noise-driven inhibition
%   alph: partial noise contribution to sustained input
%   bet:  noise effect on transient inputs
%   gon:  onset jump size. set to 0 for sutained inputs, not zero for transient or combined inputs. if nonzero, gon is computed within this program
%   goff: offset jump size. set to gon unless value provided
%   xinitial:  initial value. 0 if not provided

% outputs:
%   t: time values
%   x: firing rates
%   I: [tone noise] inputs

% example usage for Model 1 with masking and sustained inputs
%     [t,X,~] = continuityModel('sustained',1,1.5,1,3.6,5.9,1.124,.168,0,0,0,0);
%     plot(t,X)
    
% example usage for Model 2 with interrupted tones (no continuity) and transient inputs
%     [t,X,~] = continuityModel('transient',2,3,2.8,5.2,10.5,0,0,.66,1,5.19,0);
%     plot(t,X)

% example usage for Model 3 with interrupted tones and continuous responses and combined inputs
%     [t,X,~] = continuityModel('combined',2,2,1,9.5,12.7,1.124,.5,0.05,1,0.8,0);
%     plot(t,X)


%%%% Firing rate model %%%%
%%%% tau x' = -x + f(aE*x + Isustain + Ionset - Ioffset - Inoise)
%%%% where Isustain = IT + alph*IN
%%%%       Ionset   = gOn*(IT-bet*IN)*s [s exponentially decaying transient input variable]
%%%%       Ioffset  = gOff*(IT-bet*IN)*s [s exponentially decaying transient input variable]
%%%%       Inoise   = aI*IN*(1-x)

% max values
ITmax   = 5;    % maximum tone strength to allow
INmax   = 10;   % maximum noise strength to allow

%%%%% new calculation from OnePop_v3 to set tone threshold to 1 %%%%%
% m (for sustained inputs): 
k  = 1;        % steepness of nonlinearity, set this to 1

% static nonlinearity
f     = @(x) 1./(1+exp(-(x-m)/k)); % sigmoidal function

% time constants
tau     = 10; % firing rate dynamics
tauEdge = 10; % edge response expontential decay

% get jump size for transient inputs so threshold is 1 at IT=1, IN=0
if gon ~= 0

    I0 =0; switch(inputType); case('combined'); I0 = 1; case('transient0'); I0=1; end % adjust equilibrium to include sustained inputs
    [xMin,~] = fminbnd(@(x) (-x + f(aE*x+I0)), 0,1); 
    [xMax,~] = fminbnd(@(x) -(-x + f(aE*x+I0)), 0,1); 
    xR = fzero( @(x) -x + f(aE*x+I0), [1e-4,xMin-1e-4]); % rest equilibrium
    xS = fzero( @(x) -x + f(aE*x+I0), [xMin+1e-4,xMax-1e-4]); % saddle equilibrium

    % separatrix calculation.  for threshold input IT0=1
    IT0 = 1;  
    df = f(aE*xS+IT0)*(1-f(aE*xS+IT0) ) / k;
    J11 = (-1 + aE*df) / tau;
    J12 = IT0*df / tau; % without multiplier
    J22 = -1 / tauEdge;
    gon = ((J22-J11)/J12) * (xR-xS);
    
end

if isempty(goff); goff = gon; end

% set up sound start & stop info
toneStart = 200;  toneDuration = 1000; toneStop = toneStart + toneDuration;
if nTone==0 %noise alone
    noiseStart = toneStart; noiseStop = toneStop;
elseif nTone==1 % no gap
    noiseStart = toneStart; % noise starts with tone
    noiseStop = toneStop; % noise stops with tone
elseif nTone==2 % gap
    gapDuration = 500; 
    noiseStart = toneStop; noiseStop = noiseStart + gapDuration;
    toneStart(2) = noiseStop;
    toneStop(2) = toneStart(2) + toneDuration;
end

% simulation stop time
tStop  = toneStop(end) + 500;    % end time of simulation 

% functions useful for defining inputs to ODE
% define a functon that is 0 for t<T(1) and t>T(2) and 1 for T(1)<=t<=T(2)
heav1  = @(t) heaviside(t) + .5*(t==0); % 1 at origin
heav2  = @(t) heaviside(t) - .5*(t==0); % 0 at origin
rect   = @(t,T) heav1(t-T(1))-heav2(t-T(2)); % 1 at stop and start points
rect0  = @(t,T) heav2(t-T(1))-heav1(t-T(2)); % 0 at stop and start points

% sustained inputs
noiseSust = @(t) IN*rect(t,[noiseStart noiseStop]);
if nTone==0 % noise only
    toneSust = @(t) alph*IN*rect(t,[noiseStart noiseStop]);
elseif nTone==1  % no gap
    toneSust = @(t) IT*rect(t,[toneStart(1) toneStop(1)]) + alph*IN*rect(t,[noiseStart noiseStop]);
elseif nTone==2 % gap
    toneSust = @(t) IT*rect(t,[toneStart(1) toneStop(1)]) + alph*IN*rect0(t,[noiseStart noiseStop]) + ...
                    IT*rect(t,[toneStart(2) toneStop(2)]);
end

% transient inputs
if nTone==0 % noise only
    toneOn  = @(t) 0; 
    toneOff = @(t) 0;
elseif nTone==1 % tone + noise (masking)
    toneOn  = @(t) gon*exp(-(t-toneStart)/tauEdge).*heav1(t-toneStart)*max(IT-bet*IN*rect(toneStart,[noiseStart, noiseStop]), 0) ;
    toneOff = @(t) goff*exp(-(t-toneStop)/tauEdge).*heav1(t-toneStop)*max(IT-bet*IN*rect(toneStop,[noiseStart, noiseStop]), 0) ;
elseif nTone==2 % tone + noise in gap (continuity)
    toneOn  = @(t) gon*exp(-(t-toneStart(1))/tauEdge).*heav1(t-toneStart(1))*max(IT-bet*IN*rect(toneStart(1),[noiseStart, noiseStop]), 0) ...
                 + gon*exp(-(t-toneStart(2))/tauEdge).*heav1(t-toneStart(2))*max(IT-bet*IN*rect(toneStart(2),[noiseStart, noiseStop]), 0) ;
    toneOff = @(t) goff*exp(-(t-toneStop(1))/tauEdge).*heav1(t-toneStop(1))*max(IT-bet*IN*rect(toneStop(1),[noiseStart, noiseStop]), 0) ...
                 + goff*exp(-(t-toneStop(2))/tauEdge).*heav1(t-toneStop(2))*max(IT-bet*IN*rect(toneStop(2),[noiseStart, noiseStop]), 0);

end

switch(inputType)
    case('sustained')
        toneOn = @(t) zeros(size(t)); toneOff = @(t) zeros(size(t)); 
    case('transient')
        toneSust = @(t) zeros(size(t)); noiseSust = @(t) zeros(size(t));
    case('transient0')
        toneSust = @(t) zeros(size(t)); noiseSust = @(t) zeros(size(t));
    case('combined')
        % do nothing
end

% firing rate differential equation
diffEq = @(t,x)  (-x + f(aE*x + toneSust(t) + toneOn(t)-toneOff(t) - aI*(1-x)*noiseSust(t) ) )/tau;

% solve ODE 
startVec = 0; for i=1:length(toneStart); startVec(2*i:2*i+1) = [toneStart(i) toneStop(i)]; end
stopVec = [startVec(2:end) tStop];

if exist('xinitial'); x0 = xinitial;
else; x0 = 0; end

t=[]; x = [];
options = odeset('reltol',1e-8,'abstol',1e-8,'maxstep',1);

for i=1:length(startVec)
    [ti,xi] = ode15s( @(t,x) diffEq(t,x), [startVec(i) stopVec(i)], x0,options);
    t = [t; ti(1:end-1)]; x = [x; xi(1:end-1,:)];
    x0 = xi(end,:);
end

t = [t;ti(end)]; x = [x ; xi(end,:)]; 
I = [toneSust(t),toneOn(t),toneOff(t), noiseSust(t)];


