
% NOTE: Carrier multiplication not yet implemented in this (new) code
%       (will be added very soon)

figNo = 1;

% Eff vs. Eg for different etaExt
if(figNo==1)
    runSim = 1;
    if(runSim)
        Eg = linspace(0.4,2,100);
        d100 = ShockleyQueisser('Eg',Eg,'lossType','etaExt','lossVal',1);
        d90 = ShockleyQueisser('Eg',Eg,'lossType','etaExt','lossVal',0.9);
        d80 = ShockleyQueisser('Eg',Eg,'lossType','etaExt','lossVal',0.8);
        Egs = linspace(1,1.15,20);
        d100s = ShockleyQueisser('Eg',Egs,'lossType','etaExt','lossVal',1);
        d90s = ShockleyQueisser('Eg',Egs,'lossType','etaExt','lossVal',0.9);
        d80s = ShockleyQueisser('Eg',Egs,'lossType','etaExt','lossVal',0.8);
    end
    
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
elseif(figNo==2)
    runSim = 1;
    if(runSim)
        Eg = linspace(0.8,2,100);
        d100 = ShockleyQueisser('Eg',Eg,'material',@(E,Eg)StepAlpha(E,Eg,5),'lossType','etaInt','lossVal',1);
        d90 = ShockleyQueisser('Eg',Eg,'lossType','etaInt','lossVal',0.9);
        d80 = ShockleyQueisser('Eg',Eg,'lossType','etaInt','lossVal',0.8);
        d70 = ShockleyQueisser('Eg',Eg,'lossType','etaInt','lossVal',0.7);
        d1 = ShockleyQueisser('Eg',Eg,'lossType','etaInt','lossVal',0.01);
        fprintf('GaAs eff: %g \n Si eff: %g\n',dgaas.eff,dsi.eff);
    end
    
    figure;
    plot(Eg, d100.eff, 'k', Eg, d90.eff, 'b', ...
         Eg, d80.eff, 'r', Eg, d70.eff, 'g', Eg, d1.eff)
%      , ...         1.4, dgaas.eff, 'kx', 1.1, dsi.eff, 'kx');
    xlabel('Bandgap (eV)')
    ylabel('Conversion Efficiency')
    axis([0.4 2 0 35])
    legend({'\eta_{int}=100%','\eta_{int}=90%','\eta_{int}=80%','\eta_{int}=1%'});
elseif(figNo==3)
    rMin = 0.001;
    rMax = 0.999; % don't set to 1 - run into divide by zeros, etc.
    numR = 10;
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
end