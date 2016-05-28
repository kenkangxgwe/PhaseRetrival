function tarPhi = RefPhaseRetrival2(srcImg, srcZ, tarImg, tarZ, pixelsize, wavelength, FilePath)

gpuDevice(1);
resize = 6;
k = 2 * pi / wavelength;  %Wave Vector
NAs = 0.35;
resPath = [FilePath 'Result\'];

if(exist([resPath 'currentphase.mat'], 'file') )
    load([resPath 'currentphase.mat'], 'curPhi');
    srcPhi = curPhiSmoothed;
else
    error('currentphase.mat does not exist');
end

gpuSrcImg = sqrt(complex(imresize(gpuArray(srcImg), resize) ) );
[rowSize, colSize] = size(gpuSrcImg);
mesh = pixelsize / resize;

% extend the image on x direction
eySize = ceil(rowSize * 1.1);
exSize = ceil(colSize * 1.2);
kmeshx = 2 * pi / (mesh * exSize);
kmeshy = 2 * pi / (mesh * eySize);
kx = ( (1 : exSize) - exSize / 2 - 1) * kmeshx; 
ky = ( (1 : eySize) - eySize / 2 - 1) * kmeshy;
kx = fftshift(kx); ky = fftshift(ky);
[kkx, kky] = meshgrid(kx, ky);
kkx = gpuArray(kkx); kky = gpuArray(kky);
kwindow = exp(- (kkx .^ 2+kky .^ 2) / k^2 / NAs^2);
kwindow = kwindow > (max(kwindow(:) ) / 2.71828);
tarWavex = ifft2(kwindow .* exp(1i * sqrt(k^2 - kkx .^ 2 - kky .^ 2) * (tarZ - srcZ) ) .* fft2(gpuSrcImg .* exp(1i * srcPhi), eySize, exSize) );
tarWavex = tarWavex(1: rowSize , 1 : colSize);
clear kkx kky;

% extend the image on y direction
eySize = ceil(rowSize * 1.7);
exSize = ceil(colSize * 1);
mesh = pixelsize / resize;
kmeshx = 2 * pi / (mesh * exSize);
kmeshy = 2 * pi / (mesh * eySize);
kx = ( (1 : exSize) - exSize / 2 - 1) * kmeshx; 
ky = ( (1 : eySize) - eySize / 2 - 1) * kmeshy;
kx = fftshift(kx); ky = fftshift(ky);
[kkx, kky] = meshgrid(kx, ky);
kkx = gpuArray(kkx); kky = gpuArray(kky);
kwindow = exp(- (kkx .^ 2+kky .^ 2) / k^2 / NAs^2);
kwindow = kwindow > (max(kwindow(:) ) / 2.71828);
tarWavey = ifft2(kwindow .* exp(1i * sqrt(k^2 - kkx .^ 2 - kky .^ 2) * (tarZ - srcZ) ) .* fft2(gpuSrcImg .* exp(1i * srcPhi), eySize, exSize) );
tarWavey = tarWavey(1: rowSize , 1 : colSize);
clear kkx kky;

% Combine two different Images
tarWave = [tarWavey(1: rowSize, 1 : 7200), tarWavex(1: rowSize, 7201 : colSize)];
tarPhi = angle(tarWave);
tarWave = double(gather(imresize(abs(tarWave).^2, 1 / resize) ) );
figure; 
subplot(1,2,1); imagesc(abs(tarWave) ); title('propagated image'); 
subplot(1,2,2); imagesc(abs(tarImg) ); title('original image');
drawnow; 

% Calculate the residue
a = tarWave(100 : 1000, 200 : 1100); a = a(:);
b = tarImg(100 : 1000, 200 : 1100); b = b(:);

residue = 1 - (dot(a - mean(a), b - mean(b) ) / (size([a; b], 2) * std(a) * std(b) );
