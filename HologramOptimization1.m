function coef = HologramOptimization1(holI, refI, refZ, holBg, refBg, pixelsize, wavelength, FilePath)

gpuDevice(1);
resize = 6;
k = 2 * pi / wavelength;  %Wave Vector
NAs = 0.4;
crossZ = -9470; % we get the positon of cross manually
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

crossHolE = ifft2(kwindow .* exp(1i * sqrt(k^2 - kkx .^ 2 - kky .^ 2) * (crossZ - refZ) ) .* fft2(gpuHolE , eySize, exSize) );
crossRefE = ifft2(kwindow .* exp(1i * sqrt(k^2 - kkx .^ 2 - kky .^ 2) * (crossZ - refZ) ) .* fft2(gpuRefE , eySize, exSize) );
% figure; imagesc(abs(crossHolE) ); colorbar; drawnow;
% figure; imagesc(abs(crossRefE) ); colorbar; drawnow;

crossHolE =  gather(abs(crossHolE(2200 : 3600, 2400 : 3800) ) );
crossRefE =  gather(abs(crossRefE(2200 : 3600, 2400 : 3800) ) );
RoICrossHolE = abs(gradient(crossHolE) );
RoICrossRefE = abs(gradient(crossHolE) );

% Set a block at [300:480,250:750]
% block = ones(1040,1392);
% block(300:480, 250:750) = 0;
fun1 = @(x) var(crossHolE(:).^2 - x^2 * crossRefE(:).^2 );
[coef] = fminunc(fun1, 1);
coRefI = coef^2;
holIDiff = (holI - coRefI * refI);
fun1(1)
fun1(coef)
figure; imagesc(crossHolE.^2 - coRefI * crossRefE.^2); colorbar; drawnow;

fun2 = @(x) rms(holIDiff(:) - x(1)^2 * holBg(:) + x(2)^2 * refBg(:));
[coef] = fminunc(fun2, [1,1]);
coHolBg = coef(1)^2;
coRefBg = coef(2)^2;
holIDiff = (holIDiff - coHolBg * holBg - coRefBg * refBg);
fun2([1,1])
fun2(coef)

figure; imagesc(abs(holIDiff) ); colorbar; drawnow;
end