%% phaseRetrivalCaller
% A function caller, which no only calls functions step by step, but can also
% jump to a specific function directly if the program has finished the former
% functions before.

wavelength = 0.7667;
pixelsize = 6.5;
step = 4;
path = 'G:\新建文件夹\Experiment\7Apr2016\';

refZ = - 0.0505 * 25400; % Distance in mm

%% RefPhaseRetrival1
% A iterative algorithm to get the phase infomation of image a1

if(bitand(step, 1) )
    %From top to bottom, the distance decreases and the light spot shrinks
    a1 = openxlsFigures(strcat(path,'WavefrontMeasurement\'), 'waveB', 1, 100, 'up');
    bg1 = openxlsFigures(strcat(path,'WavefrontMeasurement\'), 'waveBbk', 1, 100, 'up');
    a1 = max(a1 - bg1, 0);
    z1 = - 0.2721 * 25400; 

    a2 = openxlsFigures(strcat(path,'WavefrontMeasurement\'), 'waveC', 1, 100, 'up');
    bg2 = openxlsFigures(strcat(path,'WavefrontMeasurement\'), 'waveCbk', 1, 100, 'up');
    a2 = max(a2 - bg2, 0);
    z2 = - 0.3500 * 25400;
    
    a3 = openxlsFigures(strcat(path,'WavefrontMeasurement\'), 'waveD', 1, 100, 'up');
    bg3 = openxlsFigures(strcat(path,'WavefrontMeasurement\'), 'waveDbk', 1, 100, 'up');
    a3 = max(a3 - bg3, 0);
    z3 = - 0.5000 * 25400;
    
    a4 = openxlsFigures(strcat(path,'WavefrontMeasurement\'), 'waveE', 1, 100, 'up');
    bg4 = openxlsFigures(strcat(path,'WavefrontMeasurement\'), 'waveEbk', 1, 100, 'up');
    a4 = max(a4 - bg4, 0);
    z4 = - 0.6080 * 25400;

    RefPhaseRetrival1(a4, a3, a2, a1, z4, z3, z2, z1, pixelsize, wavelength, path);
end

%% RefPhaseRetrival2
% Using angular spectrum method to get the phase infomation of reference
% image via image a1

if(bitand(step, 2) )
    refimg = openxlsFigures(strcat(path,'WavefrontMeasurement\'), 'waveA', 1, 100, 'up');
    refbg = openxlsFigures(strcat(path,'WavefrontMeasurement\'), 'waveAbk', 1, 100, 'up');
    refimg = max(refimg - refbg, 0);
    a1 = openxlsFigures(strcat(path,'WavefrontMeasurement\'), 'waveB', 1, 100, 'up');
    bg1 = openxlsFigures(strcat(path,'WavefrontMeasurement\'), 'waveBbk', 1, 100, 'up');
    a1 = max(a1 - bg1, 0);
    z1 = - 0.2500 * 25400; 

    [refPhi] = RefPhaseRetrival2(a1, z1, refimg, refZ, pixelsize, wavelength, path);
    save([path 'Result\ReferencePhase.mat'],'refPhi');
end

%% SignalPhaseRetrival
% Get the signal holograph via reference intensity, phase and signal
% intensity

if(bitand(step, 4) )
    atomCut = openxlsFigures(strcat(path,'Dipoletrap_MOT1d0V_Exp1d0ms\'), 'Dipoletrap_MOT1d0V_Exp1d0ms', 200, 200, 'up');
    atomCutRef = openxlsFigures(strcat(path,'Dipoletrap_MOT1d0V_Exp1d0ms\'), 'Dipoletrap_MOT1d0V_Exp1d0ms', 200, 200, 'down');
    atomCutBk1 = openxlsFigures(strcat(path,'Dipoletrap_MOT1d0V_Exp1d0ms\'), 'BK', 2, 10, 'up');
    atomCutBk2 = openxlsFigures(strcat(path,'Dipoletrap_MOT1d0V_Exp1d0ms\'), 'BK', 2, 10, 'down');

    SignalPhaseRetrival1(atomCut, atomCutRef, refZ, atomCutBk1, atomCutBk2, pixelsize, wavelength, path);
end