
% Need to add Geometries and Materials sub-folders to path
% No error-checking in procArgs right now

% lossType:
%   - 'none' (default)
%   - 'etaInt'
%   - 'etaExt'
%   - 'abc', [alpha beta], such that Rloss = alpha*exp(beta*qV/kT), beta
%   should be in 1/2<=beta<=3/2 (note the L-dependence, surface rate)
% material: 
%   - 'StepAlpha (default)
%   - 'GaAs'
%   - 'Si'
%   - 'lossyGaAs' (if chosen, etaInt&etaExt disabled)
%   - 'lossySi'   (if chosen, etaInt&etaExt disabled)
% geo: 
%   - 'PlaneParallel' (default)
%   - 'ppNonIdealMirror' (
%       - rearN (refractive index of material below rear interface - e.g.
%       substrate, default = 1)
%   - 'Textured'
% sun: 
%   - 'am1.5' (default)
%   - 'blackbody'

%% MAIN FUNCTION
function[data] = ShockleyQueisser(varargin)

%%%%%%%%%% defaults
sun = 'am1.5';
material = @(E,Eg) StepAlpha(E,Eg,2.5);
geo = @PlaneParallel;
lossType = 'none';
lossVal = 1;
sweepLabel = 'Eg';
Eg = linspace(0.4,2,100);
L = 1e-4;
TSun = 5780; % only if 'sun'=='blackbody'
%%%%%%%%%%

%%%%%%%%%% Fundamental Constants
q = 1.602176565e-19; 
h = 6.62606957e-34; 
kT = 0.026; % technically kT/q, approx T=302K
kB = 1.3806488e-23;
c = 29979245800; % cm/s 
%%%%%%%%%%

procArgs(varargin);

% photonFlux has units mA/cm^2/eV (i.e. flux per unit energy)
if(strcmp(sun,'am1.5'))
    fileLoc = './Data'; % location of spectrum.mat (matters if actualSun==1)
    data = load([fileLoc,'/spectrum.mat']);
    photonFlux = data.photonFlux;
    E = data.E;
    clear data;
else % this should not use the same spectrum (no need) - change at some point
    fileLoc = './Data'; % location of spectrum.mat (matters if actualSun==1)
    data = load([fileLoc,'/spectrum.mat']);
    photonFlux = data.photonFlux;
    E = data.E;
    clear data;
    omegaSun = 6.85e-5; % SQ
    photonFlux = omegaSun * 2*1e3*(q^4/h^3/c^2)*E.^2.*exp(-E/(kB*TSun/q)); % extra factor of 1e3*q to get mA
end

% totCurrent = trapz(lambda, photonFlux); % 68 mA/cm^2
totPower = trapz(E, E .* photonFlux); % 100.037 mW/cm^2 if actualSun
fprintf('Total power in solar spectrum: %g  mW/cm^2 \n', totPower);

% initialize variables
if(strcmp(sweepLabel,'Eg'))
    sweepVar = Eg;
else
    sweepVar = L;
end
numVar = length(sweepVar);
Jsc = zeros(numVar,1);
Voc = Jsc;
Vop = Jsc;
Jop = Jsc;
etaExt = Jsc;
etaInt = Jsc;
avgBounce = Jsc;
aint = Jsc;

% fast enough, no need to vectorize
for i=1:numVar  
    var = sweepVar(i); % Eg or L
    
    if( ~strcmp(lossType, 'mat') ) % a little clumsy to need this
        [alphaL, nr] = material(E, var);
    else
        [alphaL, nr, alphaBeta] = material(E, var);
    end
   [a,ar] = geo(alphaL, nr);
    
    b0 = 2*1e3*(q^4/h^3/c^2)*E.^2.*exp(-E/kT); % extra factor of 1e3*q to get mA
    intRad = 4*pi*nr^2 * trapz(E, alphaL.*b0);
    
    % compute the absorption rate
    qy = 1.*(E>=0.6) + 1.*(E>=1.2);
    Jsc(i) = trapz(E, qy.*a.*photonFlux);
    
    % compute the emission rate through the front (without exp(qV/kT)
    J0 = pi * trapz(E, a.*b0);
    Jemf = J0;   % front
    if( isa(ar,'function_handle') )
        ar = integral(@(t)ar(t).*cos(t).*sin(t)*2, 0, pi/2, 'ArrayValued', true);
    end
    Jemr = pi * trapz(E, ar.*b0); % rear
    
    % compute the loss rate
    alpha = 0; beta = 1;
    if( strcmp(lossType, 'etaExt') )
        alpha = (1-lossVal)/lossVal * J0;
        beta = 1;
    elseif( strcmp(lossType, 'etaInt') )
        alpha = (1-lossVal)/lossVal * intRad;
        beta = 1;
    elseif( strcmp(lossType, 'abc') )
        alpha = lossVal(1);
        beta = lossVal(2);
    elseif( strcmp(lossType, 'mat') )
        alpha = q * 1e3 * alphaBeta(1);
        beta = alphaBeta(2);
    end
    
    b0qy = 2*1e3*(q^4/h^3/c^2)*E.^2;
    JV = @(x) (Jsc(i) - pi * trapz(E, ar.*b0qy./(exp(E/kT-qy.*x))) ));
    dPdx = @(x) (Jsc(i) - 
%     JV = @(x) ( Jsc(i) - (Jemf+Jemr) * exp(x) - alpha*exp(beta*x) );
%     dPdx = @(x) ( Jsc(i) - (1+x) * (Jemf+Jemr) * exp(x) - (1+beta*x) * alpha*exp(beta*x) );
    xop = fzero(dPdx, 42);
    if(beta==1 || alpha==0)
        xoc = log( Jsc(i) /( Jemf + Jemr + alpha ) ); % simple formula
    else
        xoc = fzero(JV, 45);
    end
    Vop(i) = kT * xop;
    Voc(i) = kT * xoc;
    Jop(i) = JV(xop);
    
    etaExt(i) = Jemf * exp(xoc) /( Jemf*exp(xoc) + Jemr*exp(xoc) + alpha*exp(beta*xoc) );
    etaInt(i) = intRad*exp(xoc) /( intRad*exp(xoc) + alpha*exp(beta*xoc) );
    
    % where is the peak of the re-emission spectrum?
    [~,ind] = max(alphaL.*b0);
    aint(i) = 1 - a(ind) /( 4*nr^2*alphaL(ind) );
    avgBounce(i) = 1/( 4*alphaL(ind) );
end
eff = 100 * Jop .* Vop / totPower;
FF = Jop .* Vop ./ (Jsc .* Voc);

avgReabs = (1-etaInt.*aint).*etaInt.*aint ./ (1 - etaInt.*aint).^2;
avgReabsNL = aint ./ (1 - aint);
avgBounce = avgBounce .* avgReabsNL;

data = struct('Jsc',Jsc,'Voc',Voc,'Vop',Vop,...
    'Jop',Jop,'eff',eff,'FF',FF,'etaExt',etaExt,'etaInt',etaInt,...
    'avgReabs',avgReabs,'avgReabsNL',avgReabsNL,'avgBounce',avgBounce,...
    'aint',aint);
if(strcmp(sweepLabel,'Eg'))
    data.Eg = Eg;
else
    data.L = L;
end

end

%% PROCESS INPUT ARGUMENTS
function[] = procArgs(varargin)
v = varargin{:};
if( mod( length(v), 2) ~= 0 )
    error('need an even number of arguments. Exiting.');
end

% Should check this for etaInt & etaExt
%         if(val==0)
%             error('etaExt & etaInt cannot equal exactly 0.');
%         end


for i=1:2:length(v)
    var = v{i};
    val = v{i+1};
    if( strcmp(var,'sun') || strcmp(var,'Eg') || strcmp(var,'L') ...
            || strcmp(var,'TSun') || strcmp(var,'geo') ...
            || strcmp(var,'lossType') || strcmp(var,'lossVal') ...
            || strcmp(var,'rearTrans') || strcmp(var,'rearN') ) % pass straight through
        assignin('caller', var, val);
    elseif( strcmp(var,'material') )
        if( strcmp(val,'GaAs') )
            assignin('caller', var, @GaAs);
            assignin('caller', 'sweepLabel', 'L');
        elseif( strcmp(val,'Si') )
            assignin('caller', var, @Si);
            assignin('caller', 'sweepLabel', 'L');
        elseif( strcmp(val,'lossyGaAs') )
            assignin('caller', var, @GaAs);
            assignin('caller', 'sweepLabel', 'L');
            assignin('caller', 'lossType', 'mat'); 
        elseif( strcmp(val,'lossySi') )
            assignin('caller', var, @Si);
            assignin('caller', 'sweepLabel', 'L');
            assignin('caller', 'lossType', 'mat'); 
        else
            assignin('caller', var, val);
            assignin('caller', 'sweepLabel', 'Eg');
        end
    end
end

% Cannot define both etaExt and etaInt
% if lossyGaAs or lossySi, need to disable etaExt & etaInt
% 0<=rearTrans<=1
end
