function [coRefI, coHolBg, coRefBg] = AtomHologramOptimization(holI, refI, holBg, refBg, refZ, atomZ, params, FilePath)

gpuDevice(1);
wavelength =  params.wavelength;
pixelsize = params.pixelsize;
NAs = params.NAs;
resize = 6;
k = 2 * pi / wavelength;  %Wave Vector

resPath = [FilePath 'Result\'];

if(exist([resPath 'ReferencePhase.mat'], 'file') )
    load([resPath 'ReferencePhase.mat'], 'refPhi');
else
    error('ReferencePhase.mat does not exist');
end

holE = double(imresize(sqrt(holI), resize, 'nearest') );
refE = double(imresize(sqrt(refI), resize, 'nearest') ); 

[rowSize, colSize] = size(holE);
eySize = ceil(rowSize * 1);
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

gpuHolE = gpuArray(holE .* exp(1i * refPhi) );
gpuRefE = gpuArray(refE .* exp(1i * refPhi) );

atomHolE = ifft2(kwindow .* exp(1i * sqrt(k^2 - kkx .^ 2 - kky .^ 2) * (atomZ - refZ) ) .* fft2(gpuHolE , eySize, exSize) );
atomRefE = ifft2(kwindow .* exp(1i * sqrt(k^2 - kkx .^ 2 - kky .^ 2) * (atomZ - refZ) ) .* fft2(gpuRefE , eySize, exSize) );

% Craft the area to enlarge the atom point
atomHolE =  gather(abs(atomHolE(3030 : 3170, 4020 : 4120) ) );
atomRefE =  gather(abs(atomRefE(3030 : 3170, 4020 : 4120) ) );

% Optimize by minimize the variance of the difference between 
% the signal and the reference image in the range of interest
fun1 = @(x) var(atomHolE(:).^2 - x^2 * atomRefE(:).^2 );
[coef] = fminunc(fun1, 1);
coRefI = coef^2;
holIDiff = (holI - coRefI * refI);
fun1(1)
fun1(coef)
subplot(2,2,3); title(['atomHolE.^2 - coRefI * atomRefE.^2']); imagesc(atomHolE.^2 - coRefI * atomRefE.^2); colorbar;

% Optimize by minimize the root mean square of difference between
% the images and their backgrouds
fun2 = @(x) rms(holIDiff(:) - x(1)^2 * holBg(:) + x(2)^2 * refBg(:));
[coef] = fminunc(fun2, [1,1]);
coHolBg = coef(1)^2;
coRefBg = coef(2)^2;
holIDiff = (holIDiff - coHolBg * holBg - coRefBg * refBg);
fun2([1,1])
fun2(coef)

figure; title(['holIDiff']); imagesc(abs(holIDiff) ); colorbar;
end
