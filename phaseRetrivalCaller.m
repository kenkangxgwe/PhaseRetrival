%% phaseRetrivalCaller
% A function caller, which no only calls functions step by step, but can also
% jump to a specific function directly if the program has finished the former
% functions before.

wavelength = 0.7667;
pixelsize = 6.5;
step = 6;
path = 'G:\新建文件夹\Experiment\7Apr2016\';

refZ = - 0.0505 * 25400; % Distance in mm
z1 = - 0.2721 * 25400; 
z2 = - 0.3500 * 25400;
z3 = - 0.5000 * 25400;
focZ = - 0.6080 * 25400;

%% RefPhaseRetrival1
% A iterative algorithm to get the phase infomation of image a1

if(bitand(step, 1) )
    %From top to bottom, the distance decreases and the light spot shrinks
    a1 = openxlsFigures(strcat(path,'WavefrontMeasurement\'), 'waveB', 'up', 1, 100);
    bg1 = openxlsFigures(strcat(path,'WavefrontMeasurement\'), 'waveBbk', 'up', 1, 100);
    a1 = max(a1 - bg1, 0);

    a2 = openxlsFigures(strcat(path,'WavefrontMeasurement\'), 'waveC', 'up', 1, 100);
    bg2 = openxlsFigures(strcat(path,'WavefrontMeasurement\'), 'waveCbk', 'up', 1, 100);
    a2 = max(a2 - bg2, 0);
    
    a3 = openxlsFigures(strcat(path,'WavefrontMeasurement\'), 'waveD', 'up', 1, 100);
    bg3 = openxlsFigures(strcat(path,'WavefrontMeasurement\'), 'waveDbk', 'up', 1, 100);
    a3 = max(a3 - bg3, 0);
    
    a4 = openxlsFigures(strcat(path,'WavefrontMeasurement\'), 'waveE', 'up', 1, 100);
    bg4 = openxlsFigures(strcat(path,'WavefrontMeasurement\'), 'waveEbk', 'up', 1, 100);
    a4 = max(a4 - bg4, 0);

    RefPhaseRetrival1(a4, a3, a2, a1, z4, z3, z2, z1, pixelsize, wavelength, path);
end

%% RefPhaseRetrival2
% Using angular spectrum method to get the phase infomation of reference
% image via image a1

if(bitand(step, 2) )
    refimg = openxlsFigures(strcat(path,'WavefrontMeasurement\'), 'waveA', 'up', 1, 100);
    refbg = openxlsFigures(strcat(path,'WavefrontMeasurement\'), 'waveAbk', 'up', 1, 100);
    refimg = max(refimg - refbg, 0);
    a1 = openxlsFigures(strcat(path,'WavefrontMeasurement\'), 'waveB', 'up', 1, 100);
    bg1 = openxlsFigures(strcat(path,'WavefrontMeasurement\'), 'waveBbk', 'up', 1, 100);
    a1 = max(a1 - bg1, 0);

    [refPhi] = RefPhaseRetrival2(a1, z1, refimg, refZ, pixelsize, wavelength, path);
    save([path 'Result\ReferencePhase.mat'],'refPhi');
end

%% SignalPhaseRetrival
% Get the signal holograph via reference intensity, phase and signal intensity

if(bitand(step, 4) )
    atomCut = openxlsFigures(strcat(path,'OpticalLattice_MOT1d0V_Exp1d0ms\'), 'OpticalLattice_MOT1d0V_Exp1d0ms200', 'up');
    atomCutRef = openxlsFigures(strcat(path,'OpticalLattice_MOT1d0V_Exp1d0ms\'), 'OpticalLattice_MOT1d0V_Exp1d0ms200', 'down');
    atomCutBk1 = openxlsFigures(strcat(path,'OpticalLattice_MOT1d0V_Exp1d0ms\'), 'BK', 'up', 2, 10);
    atomCutBk2 = openxlsFigures(strcat(path,'OpticalLattice_MOT1d0V_Exp1d0ms\'), 'BK', 'down', 2, 10);
    [coRef, coHolbg, coRefbg] = HologramOptimization1(atomCut, atomCutRef, refZ, atomCutBk1, atomCutBk2, pixelsize, wavelength, path);
%     coRef = 1; coHolbg = 1; coRefbg = 0;
    SignalImageRetrival1(atomCut, coRef * atomCutRef, refZ, focZ, coHolbg * atomCutBk1, coRefbg * atomCutBk2, pixelsize, wavelength, path);
end