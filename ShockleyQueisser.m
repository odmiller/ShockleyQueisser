
% usage: 
%	d = ShockleyQueisser(varargin)
%	e.g. d = ShockleyQueisser('Eg',linspace(0.5,2,100),'material',@(E,Eg)StepAlpha(E,Eg,2.5));
%   look at genFigs.m for a number of examples
%
% Options:
% specify Eg (bandgap) OR L (thickness), which becomes the 'sweep variable'
% specify material as a function @(E,sweepvar)
%	e.g. stepAlpha with Eg=1.5 and sweepVar=1: @(E,L) stepAlpha(E,1.5,alpha*L);
%   or use a predefined material
%     'GaAs'      (lossless by default, can add etaInt or etaExt)
%     'Si'		  (lossless by default, can add etaInt or etaExt)
%     'lossyGaAs' (real material losses, etaInt&etaExt disabled)
%     'lossySi'   (real material losses, etaInt&etaExt disabled)
%	  @StepAlpha  (see above)
%   can also define any custom function
% specify loss by setting one of: etaInt, etaExt, abc 
%	default = lossless
%	etaExt = prob. of absorbed photon eventually escaping (i.e. LED	extraction efficiency) at open-circuit
%	etaInt = internal prob. of radiative recombination (vs. non-radiative), for a single emission event
%	for abc: ShockleyQueisser(...,'abc',[alpha beta],...)
%	    RLoss = alpha*L*exp(beta*qV/kT)
%		1/2 <= beta <= 3/2 (usually)
%		meant to make it easy to add 'ABC' coefficients
%		if sweepVar==Eg, then need to specify ('abc',[alpha*L beta])
% specify geometry, as a function
%   @PlaneParallel (default)
%   @ppNonIdealMirror
%       rearN (refractive index of material below rear interface - e.g.
%       substrate)
%   @Textured
%	or custom...
% specify sun 
%   'am1.5' (default)
%   'blackbody'
% specify concentration (defaults = 1)
%	'conc' = typical concentration through increase in short-circuit
%		current density
%	'emAngConc' = concentration through restriction of front emission angle
%		e.g. cf. "Kosten et. al., Light: Science and Applications 2, e45"
%	both values are restricted to range 1<=val<=46000
% specify carrier multiplication, 'cm', with quantum yield qy
%	cf. Hanna&Nozik, J. Appl. Phys. 100, 074510 (2006)
%	if sweepVar==Eg, expect qy(E,Eg), else qy(E)

%% MAIN FUNCTION
function[data] = ShockleyQueisser(varargin)

%%%%%%%%%% defaults
sun = 'am1.5';
material = @(E,Eg) StepAlpha(E,Eg,2.5);
geo = @PlaneParallel;
lossType = 'none';
lossVal = 1;
sweepVar = 'Eg';
sweepVal = linspace(0.4,2,100);
conc = 1;
emAngConc = 1;
qyFunc = @(E,Eg) 1;
TSun = 5780; % only used if 'sun'=='blackbody'
numE = 1000; % only used if 'sun'=='blackbody'
%%%%%%%%%%

%%%%%%%%%% Fundamental Constants
q = 1.602176565e-19; 
h = 6.62606957e-34; 
kT = 0.026; % technically kT/q, approx T=302K
kB = 1.3806488e-23;
c = 29979245800; % cm/s 
%%%%%%%%%%

procArgs(varargin);

% photonFlux has units mA/cm^2/eV (i.e. current density per unit energy)
if(strcmp(sun,'am1.5'))
    data = load('spectrum.mat'); % in ./Data
    photonFlux = data.photonFlux;
    E = data.E;
    clear data;
else 
	E = linspace(0.1, 5, numE);
    omegaSun = 6.85e-5; % SQ
    photonFlux = omegaSun * 2*1e3*(q^4/h^3/c^2)*E.^2.*exp(-E/(kB*TSun/q)); % extra factor of 1e3*q to get mA
end

% totCurrent = trapz(E, photonFlux) % 69 mA/cm^2, un-comment to verify
totPower = trapz(E, conc * E .* photonFlux); % 100.037 mW/cm^2 if 'am1.5'
fprintf('Total power in solar spectrum: %g  mW/cm^2 \n', totPower);

% initialize variables
numVar = length(sweepVal);
Jsc = zeros(numVar,1);
Voc = Jsc;
Vop = Jsc;
Jop = Jsc;
etaExt = Jsc;
etaInt = Jsc;
avgBounce = Jsc;
aint = Jsc;

% fast enough, no need to vectorize (fzero not vectorizable anyhow, I don't
%     think)
maxError1 = 0;
maxError2 = 0;
for i=1:numVar  
    val = sweepVal(i); % Eg or L
    
	[alphaL, nr] = material(E, val); % nr = refractive index
	[a, ar] = geo(alphaL, nr); % a (ar) = front (rear) absorptivity
	if( strcmp(sweepVar,'Eg') )
		qy = qyFunc(E,val);
	else
		qy = qyFunc(E);
	end
	
    b0 = 2*1e3*(q^4/h^3/c^2)*E.^2.*exp(-E/kT); % extra factor of 1e3*q to get mA
    intRad = 4*pi*nr^2 * trapz(E, qy.*alphaL.*b0);
    
    % compute the absorption rate
    Jsc(i) = conc * trapz(E, qy.*a.*photonFlux);
    
    % compute the emission rate through the front/rear - without exp(qV/kT)
    J0 = 1/emAngConc * pi * trapz(E, qy.*a.*b0);
    Jemf = J0;   % front
    if( isa(ar,'function_handle') ) % integrate over rear angle, if necessary (e.g. nonIdealMirror)
        ar = integral(@(t) ar(t).*cos(t).*sin(t)*2, 0, pi/2, 'ArrayValued', true);
    end
    Jemr = pi * trapz(E, qy.*ar.*b0); % rear
    
    % compute the internal loss rate
    alpha = 0; beta = 1;
	if( strcmp(lossType, 'etaExt') )
        alpha = (1-lossVal)/lossVal * J0;
        beta = 1;
    elseif( strcmp(lossType, 'etaInt') )
        alpha = (1-lossVal)/lossVal * intRad;
        beta = 1;
    elseif( strcmp(lossType, 'abc') )
        alpha = q * 1e3 * lossVal(1);
        beta = lossVal(2);
		if( strcmp(sweepVar,'L') )
			alpha = alpha * val;
		end
    elseif( strcmp(lossType, 'mat') )
        [~, ~, alphaBeta] = material(E, val);
        alpha = q * 1e3 * alphaBeta(1);
        beta = alphaBeta(2);
		if( strcmp(sweepVar,'L') )
			alpha = alpha * val;
		end
	end
    
	% J-V equation
    JV = @(x) ( Jsc(i) - (Jemf+Jemr) * exp(x) - alpha*exp(beta*x) );
    
	% maximize the power P=J*V, set dP/dx = 0
	dPdx = @(x) ( Jsc(i) - (1+x) * (Jemf+Jemr) * exp(x) - (1+beta*x) * alpha*exp(beta*x) );
    xop = fzero(dPdx, 42);
    if(beta==1 || alpha==0)
        xoc = log( Jsc(i) /( Jemf + Jemr + alpha ) ); % simple formula
    else
        xoc = fzero(JV, 45);
    end
    Vop(i) = kT * xop;
    Voc(i) = kT * xoc;
    Jop(i) = JV(xop);

	% etaExt and etaInt at open-circuit
	etaExt(i) = Jemf * exp(xoc) /( Jemf*exp(xoc) + Jemr*exp(xoc) + alpha*exp(beta*xoc) );
    etaInt(i) = intRad*exp(xoc) /( intRad*exp(xoc) + alpha*exp(beta*xoc) );
    
    % where is the peak of the re-emission spectrum?
	%	not correct for qy(E)~=1
    [~,ind] = max(alphaL.*b0);
    aint(i) = 1 - a(ind) /( 4*nr^2*alphaL(ind) );
    avgBounce(i) = 1/( 4*alphaL(ind) );
	
	% integral error estimates
	int1 = trapz(E,a.*E.^2.*exp(xop-E/kT));
	int2 = trapz(E,a.*E.^2.*exp(2*(xop-E/kT))); % 1st correction in Taylor expansion
	int3 = trapz(E, a.*qy.*E.^2./(exp(E/kT-qy*xop)));
	int4 = trapz(E, a.*qy.*E.^2./(exp(E/kT-xop)));
	if(int2/int1 > maxError1)
		maxError1 = int2/int1;
	end
	if((int4-int3)/int4 > maxError2)
		maxError2 = (int4-int3)/int4;
	end
end
fprintf('Rough estimates of max relative approximation errors \n');
fprintf('    Boltzmann approximation: %g \n', maxError1);
fprintf('    QY=1 in integral approx: %g \n', maxError2);

eff = 100 * Jop .* Vop / totPower;
FF = Jop .* Vop ./ (Jsc .* Voc);

avgReabs = (1-etaInt.*aint).*etaInt.*aint ./ (1 - etaInt.*aint).^2;
avgReabsNL = aint ./ (1 - aint);
avgBounce = avgBounce .* avgReabsNL;

data = struct('Jsc',Jsc,'Voc',Voc,'Vop',Vop,...
    'Jop',Jop,'eff',eff,'FF',FF,'etaExt',etaExt,'etaInt',etaInt,...
    'avgReabs',avgReabs,'avgReabsNL',avgReabsNL,'avgBounce',avgBounce,...
    'aint',aint,'sweepVar',sweepVar,'sweepVal',sweepVal);

end

%% PROCESS INPUT ARGUMENTS
% error checking:
%   length(varargin) == even
%   0 < etaInt, etaInt <= 1
%	1 <= conc, emAngCong <= 46000
%   numE >= 1
%   material function has exactly two inputs
%   geo function has exactly two inputs
%	can't define both Eg and L
%   can't define both etaInt and etaExt
%	can't define lossy material and etaInt/etaExt
%   if 'abc', length(val)==2, val(1)>0, 0.5 <= val(2)=beta <= 1.5
%   Note: I don't check for re-defined variables (last one wins)
function[] = procArgs(varargin)
v = varargin{:};
if( mod( length(v), 2) ~= 0 )
    error('need an even number of arguments');
end

% check for clashing args
lossyMatDef = cellfun(@(x)(strcmp(x,'lossyGaAs')+strcmp(x,'lossySi')), v);
etaIntDef = cellfun(@(x) strcmp(x,'etaInt'), v);
etaExtDef = cellfun(@(x) strcmp(x,'etaExt'), v);
if( sum(lossyMatDef)>0 && ( sum(etaIntDef)>0 || sum(etaExtDef)>0 ) )
	error('Cannot define a lossy material AND independent loss parameters');
end
if( sum(etaIntDef)>0 && sum(etaExtDef)>0 )
	error('Cannot define both etaInt and etaExt');
end

EgDef = cellfun(@(x) strcmp(x,'Eg'), v);
LDef = cellfun(@(x) strcmp(x,'L'), v);
if( sum(EgDef)>0 && sum(LDef)>0 )
	error('Cannot define both L and Eg');
end

for i=1:2:length(v)
    var = v{i};
    val = v{i+1};
    if( strcmp(var,'sun') || strcmp(var,'TSun') ...
            || strcmp(var,'rearTrans') || strcmp(var,'rearN')   ) % pass straight through
        assignin('caller', var, val);
	elseif( strcmp(var,'geo') )
		if( nargin(val)~=2 )
			error('geo function must have exactly two arguments');
		end
		assignin('caller', var, val);
	elseif( strcmp(var,'conc') || strcmp(var,'emAngConc')  ) % pass straight through
		if( val<1 || val>46000 )
			error('conc or emAngConc out of bounds [1,46000]');
		end
		assignin('caller', var, val);
	elseif( strcmp(var,'numE') )
		if(val<1)
			error('Need at least one point in solar spectrum (numE>=1)');
		end
		assignin('caller',var,val);
	elseif( strcmp(var,'cm') )
		assignin('caller', 'qyFunc', val);
	elseif( strcmp(var,'etaExt') || strcmp(var,'etaInt') )
		if( val<=0 || val>1 )
			error('etaInt or etaExt out of bounds (0,1]');
		end
		assignin('caller', 'lossType', var);
		assignin('caller', 'lossVal', val);
	elseif( strcmp(var,'abc') )
		if( length(val)~=2 )
			error('abc losses require a length-two vector argument');
		end
		if( val(1)<=0 || val(2)<0.5 || val(2)>1.5 )
			error('abc vector out of bounds: alpha>0, 0.5<=beta<=1.5');
		end
		assignin('caller', 'lossType', var);
		assignin('caller', 'lossVal', val);
	elseif( strcmp(var,'Eg') || strcmp(var,'L') )
		assignin('caller', 'sweepVar', var);
		assignin('caller', 'sweepVal', val);
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
		if( nargin(evalin('caller',var)) ~= 2 )
			error('Must define a material function with two input args');
		end
	else
		warning('Unknown option: %s.  Ignoring...', var);
    end
end

end

