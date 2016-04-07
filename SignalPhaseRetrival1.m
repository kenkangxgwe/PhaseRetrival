function SignalPhaseRetrival1(holI, refI, refZ, holBg, refBg, pixelsize, wavelength, FilePath)

gpuDevice(1);
resize = 6;
k = 2 * pi / wavelength;  %Wave Vector
NAs = 0.4;
resPath = [FilePath 'Result\'];

if(exist([resPath 'ReferencePhase.mat'], 'file') )
    load([resPath 'ReferencePhase.mat'], 'refPhi');
else
    error('ReferencePhase.mat does not exist');
end

% Set a block at [300:480,250:750]
block = ones(1040,1392);
block(300:480, 250:750) = 0;

% fun1 = @(x) rms( (holI (block == 1) - x^2 * refI (block == 1) ) );
% [coef] = fminunc(fun1, 1);
% coRefI = coef^2;
coRefI = 1;
holDiff = (holI - coRefI * refI);
refDiff = refI - refBg;
holDiff = holDiff .* (refDiff >= 0.01);
refDiff = max(refDiff, 0.01);

% fun2 = @(x) rms(rms(holDiff - x(1)^2 * holBg + x(2)^2 * refBg) );
% [coef] = fminunc(fun2,[1,1]);
% coHolBg = coef(1)^2;
% coRefBg = coef(2)^2;
% holDiff = holDiff - coHolBg * holBg + coRefBg * refBg;

% fun = @(x) rms(rms( (holDiff - x(1)^2 * holBg + x(2)^2 * refBg) .* (holI <= 1020) ) );
% [coef] = fminunc(fun2,[1,1]);
% coHolBg = coef(1)^2;
% coRefBg = coef(2)^2;
% holDiff = holDiff - coHolBg * holBg + coRefBg * refBg;
% holDiff = abs(holDiff);

holDiff = double(imresize(holDiff, resize, 'nearest') );
refE = double(imresize(sqrt(refDiff), resize, 'nearest') ); % E = sqrt(I)

holE = holDiff./ (refE .* exp(- 1i * refPhi) );

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

nxDisp=2000:4000;
nyDisp=2000:4000;
gpuHolE = gpuArray(holE);

figure;
% for sigZ = (- 0.0505 * 25400) : -100 : (- 0.6050 * 25400)
    sigZ = -0.573 * 25400;
    sigE = ifft2(kwindow .* exp(1i * sqrt(k^2 - kkx .^ 2 - kky .^ 2) * (sigZ - refZ) ) .* fft2(gpuHolE , eySize, exSize) );
    imagesc(abs(sigE(nyDisp,nxDisp) ) ); colorbar; title(['sigZ = ' num2str(sigZ)]);
    drawnow;
%      spaceMat(:, :, z) = abs(sigE(x,y) );
end
