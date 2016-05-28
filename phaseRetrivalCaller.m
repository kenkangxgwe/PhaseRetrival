function phaseRetrivalCaller(confFile)
% A function caller, all parameters are set by configuration .xml file.
% The function not only calls functions step by step,
% but can also jump to a specific function directly
% if the program has finished the former functions before.
allparams = xmlread(confFile);
imageparams = allparams.getSoleElement('images');
focusimage = imageparams.getSoleElement('focusimage');
medianimage1 = imageparams.getSoleElement('medianimage1');
medianimage2 = imageparams.getSoleElement('medianimage2');
medianimage3 = imageparams.getSoleElement('medianimage3');
referenceimage = imageparams.getSoleElement('referenceimage');

wavelength =  allparams.getSoleElementData('wavelength').str2double;
pixelsize = allparams.getSoleElementData('pixelsize').str2double;
step = bin2dec(allparams.getSoleElementData('step').char);
path = allparams.getSoleElementData('path').char;

% Distance in um
focZ = 0.625 * focusimage.getSoleElementData('z').str2double;
medZ1 = 0.625 * medianimage1.getSoleElementData('z').str2double; 
medZ2 = 0.625 * medianimage2.getSoleElementData('z').str2double; 
medZ3 = 0.625 * medianimage3.getSoleElementData('z').str2double; 
refZ = 0.625 * referenceimage.getSoleElementData('z').str2double; 


%% RefPhaseRetrival1
% A iterative algorithm to get the phase infomation of image a1

if(bitand(step, 1) )
    focusimagebg = imageparams.getSoleElement('focusimagebg');
    medianimagebg1 = imageparams.getSoleElement('medianimagebg1');
    medianimagebg2 = imageparams.getSoleElement('medianimagebg2');
    medianimagebg3 = imageparams.getSoleElement('medianimagebg3');

    % Loading the image from smail to big

    focImg = openxmlFigures(path, focusimage);
    focBg = openxmlFigures(path, focusimagebg);
    focImg = max(focImg - focBg, 0);

    medImg1 = openxmlFigures(path, medianimage1);
    medBg1 = openxmlFigures(path, medianimagebg1);
    medImg1 = max(medImg1 - medBg1, 0);
    
    medImg2 = openxmlFigures(path, medianimage2);
    medBg2 = openxmlFigures(path, medianimagebg2);
    medImg2 = max(medImg2 - medBg2, 0);

    medImg3 = openxmlFigures(path, medianimage3);
    medBg3 = openxmlFigures(path, medianimagebg3);
    medImg3 = max(medImg3 - medBg3, 0);

    RefPhaseRetrival1(focImg, medImg1, medImg2, medImg3, focZ, medZ1, medZ2, medZ3, pixelsize, wavelength, path);
end

%% RefPhaseRetrival2
% Using angular spectrum method to get the phase infomation of refImg via medImg3

if(bitand(step, 2) )
    medianimagebg3 = imageparams.getSoleElement('medianimagebg3');
    referenceimagebg = imageparams.getSoleElement('referenceimagebg');
    
    medImg3 = openxmlFigures(path, medianimage3);
    medBg3 = openxmlFigures(path, medianimagebg3);
    medImg3 = max(medImg3 - medBg3, 0);
    
    refImg =  openxmlFigures(path, referenceimage);
    refbg = openxmlFigures(path, referenceimagebg);
    refImg = max(refImg - refbg, 0);

    [refPhi] = RefPhaseRetrival2(medImg3, medZ3, refImg, refZ, pixelsize, wavelength, path);
    save([path 'Result\ReferencePhase.mat'],'refPhi');
end

%% SignalPhaseRetrival
% Get the signal holograph via reference intensity, phase and signal intensity

if(bitand(step, 4) )

    atomimage = imageparams.getSoleElement('atomimage');
    atomimagebg = imageparams.getSoleElement('atomimagebg');
    
    [atomImg] = openxmlFigures(path, atomimage);
    atomRefImg = atomImg(:,:,2);
    atomImg = atomImg(:,:,1);
    [atomBg] = openxmlFigures(path, atomimagebg);
    atomRefBg = atomBg(:,:,2);
    atomBg = atomBg(:,:,1);
    
    [coRef, coHolbg, coRefbg] = AtomHologramOptimization(atomImg, atomRefImg, atomBg, atomRefBg, refZ, pixelsize, wavelength, path);
    % coRef = 1; coHolbg = 1; coRefbg = 0;
    SignalImageRetrival1(atomCut, coRef * atomCutRef, coHolbg * atomCutBk1, coRefbg * atomCutBk2, refZ, focZ, pixelsize, wavelength, path);
end

end

%% Functions dealing with xml file

function child = getSoleElement(parent, childname)
child = parent.getElementsByTagName(childname).item(0);
end

function  childData= getSoleElementData(parent, childname)
childData = parent.getElementsByTagName(childname).item(0).getTextContent;
end

function figure = openxmlFigures(path, image)
    FilePath = image.getSoleElementData('filepath').char;
    FileName = image.getSoleElementData('filename').char;
    if(isempty(image.getSoleElement('pos') ) )
        if(isempty(image.getSoleElement('startindex') ) )
            figure1 = openxlsFigures([path FilePath], FileName, 'up');
            figure2 = openxlsFigures([path FilePath], FileName, 'down');
        else
            startindex = image.getSoleElementData('startindex').str2num;
            endindex = image.getSoleElementData('endindex').str2num;
            figure1 = openxlsFigures([path FilePath], FileName, 'up', startindex, endindex);
            figure2 = openxlsFigures([path FilePath], FileName, 'down', startindex, endindex);
        end
        figure(:,:,1) = figure1;
        figure(:,:,2) = figure2;
    else
        pos = image.getSoleElementData('pos').char;
        startindex = image.getSoleElementData('startindex').str2num;
        endindex = image.getSoleElementData('endindex').str2num;
        figure = openxlsFigures([path FilePath], FileName, pos, startindex, endindex);
    end
end