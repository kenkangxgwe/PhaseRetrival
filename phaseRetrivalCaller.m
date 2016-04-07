%% phaseRetrivalCaller
% A function caller, which no only calls functions step by step, but can also
% jump to a specific function directly if the program has finished the former
% functions before.

wavelength = 0.7667;
pixelsize = 6.5;
step = 2;
path = 'G:\新建文件夹\Experiment\20Mar2016\';

% bg1 = openxlsFigures(path, 'BK', 1, 10, 'up');
% bg2 = openxlsFigures(path, 'BK', 1, 10, 'down');
refZ = - 0.0505 * 25400; % Distance in mm

%% RefPhaseRetrival1
% A iterative algorithm to get the phase infomation of image a1

if(bitand(step, 1) )
    %From top to bottom, the distance decreases and the light spot shrinks
    a1 = openxlsFigures(path, 'waveB', 1, 100, 'up');
    bg1 = openxlsFigures(path, 'waveBbk', 1, 100, 'up');
    a1 = max(a1 - bg1, 0);
    z1 = - 0.2500 * 25400; 

    a2 = openxlsFigures(path, 'waveC', 2, 100, 'up');
    bg2 = openxlsFigures(path, 'waveCbk', 23, 100, 'up');
    a2 = max(a2 - bg2, 0);
    z2 = - 0.3500 * 25400;
    
    a3 = openxlsFigures(path, 'waveD', 1, 100, 'up');
    bg3 = openxlsFigures(path, 'waveDbk', 1, 100, 'up');
    a3 = max(a3 - bg3, 0);
    z3 = - 0.5000 * 25400;
    
    a4 = openxlsFigures(path, 'waveE', 1, 100, 'up');
    bg4 = openxlsFigures(path, 'waveEbk', 2, 89, 'up');
    a4 = max(a4 - bg4, 0);
    z4 = - 0.6050 * 25400;

    RefPhaseRetrival1(a4, a3, a2, a1, z4, z3, z2, z1, pixelsize, wavelength, path);
end

%% RefPhaseRetrival2
% Using angular spectrum method to get the phase infomation of reference
% image via image a1

if(bitand(step, 2) )
    refimg = openxlsFigures(path, 'waveA', 1, 100, 'up');
    refbg = openxlsFigures(path, 'waveAbk', 1, 100, 'up');
    refimg = max(refimg - refbg, 0);
    a1 = openxlsFigures(path, 'waveB', 1, 100, 'up');
    bg1 = openxlsFigures(path, 'waveBbk', 1, 100, 'up');
    a1 = max(a1 - bg1, 0);
    z1 = - 0.2500 * 25400; 

    [refPhi] = RefPhaseRetrival2(a1, z1, refimg, refZ, pixelsize, wavelength, path);
    save([path 'Result\ReferencePhase.mat'],'refPhi');
end

%% SignalPhaseRetrival
% Get the signal holograph via reference intensity, phase and signal
% intensity

if(bitand(step, 4) )
%     atomCut = openxlsFigures(path, 'p10n10acut', 8, 15, 'up');
%     atomCutRef = openxlsFigures(path, 'p10n10acut', 8, 15, 'down');
    atomCut = openxlsFigures(path, 'fscancut', 1, 1, 'up');
    atomCutRef = openxlsFigures(path, 'fscancut', 1, 1, 'down');
    atomCutBk1 = openxlsFigures(path, 'fscannolabk', 1, 6, 'up');
    atomCutBk2 = openxlsFigures(path, 'fscannolabk', 1, 6, 'down');
%     atomNoCut = openxlsFigures(path, 'p10n10acut', 14, 14, 'up');
%     at    omNoCutRef = openxlsFigures(path, 'p10n10acut', 14, 14, 'down');

%     atomNoCut = max(atomNoCut - bg1, 10);
%     atomNoCutRef = max(atomNoCutRef - bg1, 10);
%     atomCut = max(atomCut - bg1, 10);
%     atomCutRef = max(atomCutRef - bg1, 10);

    SignalPhaseRetrival1(atomCut, atomCutRef, refZ, atomCutBk1, atomCutBk2, pixelsize, wavelength, path);
end