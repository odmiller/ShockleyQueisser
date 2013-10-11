
% NOTE: Carrier multiplication not yet implemented in this (new) code
%       (will be added very soon)

figNo = [4]; % vector of all figNo to run

% Eff vs. Eg for different etaExt
if( ismember(1,figNo) )
	Eg = linspace(0.4,2,100);
	d100 = ShockleyQueisser('Eg',Eg,'etaExt',1);
	d90 = ShockleyQueisser('Eg',Eg,'etaExt',0.9);
	d80 = ShockleyQueisser('Eg',Eg,'etaExt',0.8);
	Egs = linspace(1,1.15,20);
	d100s = ShockleyQueisser('Eg',Egs,'etaExt',1);
	d90s = ShockleyQueisser('Eg',Egs,'etaExt',0.9);
	d80s = ShockleyQueisser('Eg',Egs,'etaExt',0.8);
    
    figure;
    plot(Eg, d100.eff, Eg, d90.eff, Eg, d80.eff)
    xlabel('Bandgap (eV)')
    ylabel('Conversion Efficiency')
    axis([0.4 2 0 35])
    legend({'\eta_{ext}=1','\eta_{ext}=0.9','\eta_{ext}=0.8'});
    
    figure;
    plot(Egs, d100s.eff, Egs, d90s.eff, Egs, d80s.eff)
    xlabel('Bandgap (eV)')
    ylabel('Conversion Efficiency')
end

% Eff vs. Eg for different etaInt
if( ismember(2,figNo) )
	Eg = linspace(0.8,2,100);
	d100 = ShockleyQueisser('Eg',Eg,'material',@(E,Eg)StepAlpha(E,Eg,5),'etaInt',1);
	d90 = ShockleyQueisser('Eg',Eg,'etaInt',0.9);
	d80 = ShockleyQueisser('Eg',Eg,'etaInt',0.8);
	d70 = ShockleyQueisser('Eg',Eg,'etaInt',0.7);
	d1 = ShockleyQueisser('Eg',Eg,'etaInt',0.01);
    
    figure;
    plot(Eg, d100.eff, 'k', Eg, d90.eff, 'b', ...
         Eg, d80.eff, 'r', Eg, d70.eff, 'g', Eg, d1.eff);
    xlabel('Bandgap (eV)')
    ylabel('Conversion Efficiency')
    axis([0.4 2 0 35])
    legend({'\eta_{int}=100%','\eta_{int}=90%','\eta_{int}=80%','\eta_{int}=1%'});
end

% Eff vs. rear mirror reflectivity for 3um thickness, GaAs
if( ismember(3,figNo) )
    rMin = 0.001;
    rMax = 0.999; % don't set to 1 - run into divide by zeros, etc.
    numR = 50;
    r = linspace(rMin,rMax,numR).';
    eff = zeros(size(r));
    voc = eff;
    jsc = eff;
	for i=1:numR
        dgaas = ShockleyQueisser('material','lossyGaAs',...
            'L',3e-4,'geo',@(al,nr)ppNonIdealMirror(al,nr,1-r(i),nr));
        eff(i) = dgaas.eff;
        voc(i) = dgaas.Voc;
        jsc(i) = dgaas.Jsc;
	end
	figure;
	plot(r,eff);
	xlabel('Reflectivity');
	ylabel('Cell Efficiency (%)');
	
	figure; 
	plot(r,voc);
	xlabel('Reflectivity');
	ylabel('V_{OC} (Volts)');
	ylim([1.06 1.16]);
	
	figure;
	plot(r,jsc);
	xlabel('Reflectivity');
	ylabel('J_{SC} (mA/cm^2)');
end

% 6x carrier multiplication, step function, for different etaInt
if( ismember(4,figNo) )
	Eg = linspace(0.4,2,100);
	qy = @(E,Eg) floor(E/Eg);
	d100 = ShockleyQueisser('Eg',Eg,'material',@(E,Eg) StepAlpha(E,Eg,5),'cm',qy);
	d90 = ShockleyQueisser('Eg',Eg,'material',@(E,Eg) StepAlpha(E,Eg,5),'cm',qy, 'etaInt',0.9);
	d80 = ShockleyQueisser('Eg',Eg,'material',@(E,Eg) StepAlpha(E,Eg,5),'cm',qy, 'etaInt',0.8);
	figure; 
	plot(Eg, d100.eff, Eg, d90.eff, Eg, d80.eff);
	xlabel('E_G (eV)');
	ylabel('Solar Cell Efficiency (%)');
	legend({'\eta_{ext}=1','\eta_{ext}=0.9','\eta_{ext}=0.8'});
end

% example of Eff. vs. Thickness for diff. etaInt, with 100x concentration
if( ismember(5,figNo) )
	L = linspace(0.4,2,100); % arbitrary units once multiplied by alpha.  Could be microns
	mat = @(E,L) StepAlpha(E,0.7,2*L); % alpha = 2/micron, say
	d100 = ShockleyQueisser('L',L,'etaInt',1,'material',mat,'conc',100);
	d90 = ShockleyQueisser('L',L,'etaInt',0.9,'material',mat,'conc',100);
% 	d100 = ShockleyQueisser('L',L,'etaInt',1,'material',mat,'emAngConc',100);
% 	d90 = ShockleyQueisser('L',L,'etaInt',0.9,'material',mat,'emAngConc',100);
	
    figure;
    plot(L, d100.eff, L, d90.eff)
    xlabel('Thickness (\mum)')
    ylabel('Conversion Efficiency')
%     axis([0.4 2 0 35])
    legend({'\eta_{ext}=1','\eta_{ext}=0.9'});
end

% check ABC def's by verifying with Si
if( ismember(6,figNo) )
	L = logspace(1,3,100)*1e-4; % cm, NOT arbitrary units
	dsi = ShockleyQueisser('L',L,'material','Si','geo',@Textured);
	dsiL = ShockleyQueisser('L',L,'material','lossySi','geo',@Textured);
	alphaBeta = [3.88E-31*1.45e10^3 3/2];
	dabc = ShockleyQueisser('L',L,'material','Si','abc',alphaBeta,'geo',@Textured);
	
	figure; 
	semilogx(L*1e4, dsi.eff, L*1e4, dsiL.eff, L*1e4, dabc.eff);
	xlabel('Thickness (\mum)');
	ylabel('Conversion Efficiency');
	
	figure;
	semilogx(L*1e4, dsi.Jsc, L*1e4, dsiL.Jsc, L*1e4, dabc.Jsc);
	xlabel('Thickness (\mum)');
	ylabel('J_{SC}');
	
	figure;
	semilogx(L*1e4, dsi.Voc, L*1e4, dsiL.Voc, L*1e4, dabc.Voc);
	xlabel('Thickness (\mum)');
	ylabel('V_{OC}');
end

% GaAs, Planar+GoodMirror vs. Textured+GoodMirror vs. Textured+BadMirror
if( ismember(7,figNo) )
	L = logspace(-1,2,100)*1e-4; % cm, NOT arbitrary units
	dgaas1 = ShockleyQueisser('L',L,'material','lossyGaAs','geo',@PlaneParallel);
	dgaas2 = ShockleyQueisser('L',L,'material','lossyGaAs','geo',@Textured);
	
	geo3 = @(al,nr)ppNonIdealMirror(al,nr,0.999,nr); % 99% transmission, substrate underneath
	dgaas3 = ShockleyQueisser('L',L,'material','lossyGaAs','geo', geo3);
	
	figure; 
	semilogx(L*1e4, dgaas1.eff, L*1e4, dgaas2.eff, L*1e4, dgaas3.eff);
	xlabel('Thickness (\mum)');
	ylabel('Conversion Efficiency');
	
	figure;
	semilogx(L*1e4, dgaas1.Jsc, L*1e4, dgaas2.Jsc, L*1e4, dgaas3.Jsc);
	xlabel('Thickness (\mum)');
	ylabel('J_{SC}');
	
	figure;
	semilogx(L*1e4, dgaas1.Voc, L*1e4, dgaas2.Voc, L*1e4, dgaas3.Voc);
	xlabel('Thickness (\mum)');
	ylabel('V_{OC}');
end
