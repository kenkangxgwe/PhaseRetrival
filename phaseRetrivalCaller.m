%% phaseRetrivalCaller
% A function caller, which not only calls functions step by step, but can also
% jump to a specific function directly if the program has finished the former
% functions before.

wavelength = 0.7667;
pixelsize = 6.5;
step = bin2dec('001');
path = 'G:\新建文件夹\Experiment\21May2016\';

refZ = 19600 * 0.625; % A
z1 = 6200 * 0.625; % B5
z2 = -300 * 0.625; % B3
z3 = - 5300 * 0.625; %B1
focZ = - 7800 * 0.625; % E

%% RefPhaseRetrival1
% A iterative algorithm to get the phase infomation of image a1

if(bitand(step, 1) )
    
% Loading the image in the order of the size from big to smail
    a1 = openxlsFigures([path 'WavefrontMeasurement\'], 'waveB5_p6200a', 'up', 1, 100);
    bg1 = openxlsFigures([path 'WavefrontMeasurement\'], 'waveB5bk_p6200a', 'up', 1, 100);
    a1 = max(a1 - bg1, 0);

    a2 = openxlsFigures([path 'WavefrontMeasurement\'], 'waveB3_n300a', 'up', 1, 100);
    bg2 = openxlsFigures([path 'WavefrontMeasurement\'], 'waveB3bk_n300a', 'up', 1, 100);
    a2 = max(a2 - bg2, 0);
    
    a3 = openxlsFigures([path 'WavefrontMeasurement\'], 'waveB1_n5300a', 'up', 1, 100);
    bg3 = openxlsFigures([path 'WavefrontMeasurement\'], 'waveB1bk_n5300a', 'up', 1, 100);
    a3 = max(a3 - bg3, 0);

    a4 = openxlsFigures([path 'WavefrontMeasurement\'], 'waveE_n7800a', 'up', 1, 100);
    bg4 = openxlsFigures([path 'WavefrontMeasurement\'], 'waveEbk_n7800a', 'up', 1, 100);
    a4 = max(a4 - bg4, 0);

    RefPhaseRetrival1(a4, a3, a2, a1, focZ, z3, z2, z1, pixelsize, wavelength, path);
end

%% RefPhaseRetrival2
% Using angular spectrum method to get the phase infomation of reference
% image via image a1

if(bitand(step, 2) )
    refimg = openxlsFigures([path 'WavefrontMeasurement\'], 'waveA_p19600a', 'up', 1, 100);
    refbg = openxlsFigures([path 'WavefrontMeasurement\'], 'waveAbk_p19600a', 'up', 1, 100);
    refimg = max(refimg - refbg, 0);
    a1 = openxlsFigures([path 'WavefrontMeasurement\'], 'waveB5_p6200a', 'up', 1, 100);
    bg1 = openxlsFigures([path 'WavefrontMeasurement\'], 'waveB5_p6200a', 'up', 1, 100);
    a1 = max(a1 - bg1, 0);

    [refPhi] = RefPhaseRetrival2(a1, z1, refimg, refZ, pixelsize, wavelength, path);
    save([path 'Result\ReferencePhase.mat'],'refPhi');
end

%% SignalPhaseRetrival
% Get the signal holograph via reference intensity, phase and signal intensity

if(bitand(step, 4) )

     atomCut = openxlsFigures(strcat(path,'OpticalLattice_MOT1d0V_Exp0d2ms\'), 'OpticalLattice_MOT1d0V_Exp0d2ms10', 'up');
     atomCutRef = openxlsFigures(strcat(path,'OpticalLattice_MOT1d0V_Exp0d2ms\'), 'OpticalLattice_MOT1d0V_Exp0d2ms10', 'down');
     atomCutBk1 = openxlsFigures(strcat(path,'OpticalLattice_MOT1d0V_Exp0d2ms\'), 'BK', 'up', 1, 21);
     atomCutBk2 = openxlsFigures(strcat(path,'OpticalLattice_MOT1d0V_Exp0d2ms\'), 'BK', 'down', 1, 21);

    [coRef, coHolbg, coRefbg] = HologramOptimization1(atomCut, atomCutRef, refZ, atomCutBk1, atomCutBk2, pixelsize, wavelength, path);
%     coRef = 1; coHolbg = 1; coRefbg = 0;
    SignalImageRetrival1(atomCut, coRef * atomCutRef, refZ, focZ, coHolbg * atomCutBk1, coRefbg * atomCutBk2, pixelsize, wavelength, path);
end