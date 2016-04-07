function tarPhi = RefPhaseRetrival2(srcimg, srcZ, tarimg, tarZ, pixelsize, wavelength, FilePath)

gpuDevice(1);
resize = 6;
k = 2 * pi / wavelength;  %Wave Vector
NAs = 0.4;
resPath = [FilePath 'Result\'];

if(exist([resPath 'currentPhase.mat'], 'file') )
    load([resPath 'currentPhase.mat'], 'curPhi');
    srcPhi = curPhi;
else
    error('currentPhase.mat does not exist');
end

srcimg = double(imresize(srcimg, resize, 'nearest') );

% Denoising, threshold is max(srcimg(:) ) / 300
srcimg = srcimg - max(srcimg(:) ) / 300;
srcimg = sqrt(srcimg .* (srcimg>0) ); % E field = sqrt(Amplitude)

[rowSize, colSize] = size(srcimg);
eySize = ceil(rowSize * 1.4);
exSize = ceil(colSize * 1);

mesh = pixelsize / resize;
kmeshx = 2 * pi / (mesh * exSize);
kmeshy = 2 * pi / (mesh * eySize);
kx = ( (1 : exSize) - exSize / 2 - 1) * kmeshx; 
ky = ( (1 : eySize) - eySize / 2 - 1) * kmeshy;
kx = fftshift(kx);
ky = fftshift(ky);
[kkx, kky] = meshgrid(kx, ky);
kkx = gpuArray(kkx);
kky = gpuArray(kky);
kwindow = exp(- (kkx .^ 2+kky .^ 2) / k^2 / NAs^2);
kwindow = kwindow > (max(kwindow(:) ) / 2.71828);

gpuSrcimg = gpuArray(srcimg);
tarWave = ifft2(kwindow .* exp(1i * sqrt(k^2 - kkx .^ 2 - kky .^ 2) * (tarZ - srcZ) ) .* fft2(gpuSrcimg .* exp(1i * srcPhi), eySize, exSize) );
tarWave = [tarWave(1: rowSize , 1 : 7500), zeros(rowSize, colSize - 7500)];
tarWave = gather(abs(tarWave) );
tarWave = double(imresize(tarWave, 1 / resize, 'nearest') );
figure; imagesc(abs(tarWave) );
residue = mean(abs( abs( tarWave(tarimg ~= 0) ) / mean( abs( tarWave(tarimg ~= 0) ) ) - abs(tarimg(tarimg ~= 0) / mean(tarimg(tarimg ~= 0) ) ) ) ./ (abs(tarimg(tarimg ~= 0) ) / mean(tarimg(tarimg ~= 0) ) ) )
tarPhi = angle(tarWave);
