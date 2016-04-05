function RefPhaseRetrival1(image0, image1, image2, image3, z0, z1, z2, z3, pixelsize, wavelength, FilePath)
% image0 focal
% image1 small
% image2 medium
% image3 big, but not being clipped by the camera
% after we get the phase for image3, we propagate E3 to the plane where we
% are going to do the experiment.

gpuDevice(1);
resize = 6;
k = 2 * pi / wavelength;  %Wave Vector
NAs = 0.4;
resPath = [FilePath '\Result\'];

image0 = double(imresize(image0, resize, 'nearest') );
image1 = double(imresize(image1, resize, 'nearest') );
image2 = double(imresize(image2, resize, 'nearest') );
image3 = double(imresize(image3, resize, 'nearest') );

%Denoising, threshold = max(nimage(:) )/300
image1 = image1 - max(image1(:) ) / 300;
image1 = sqrt(image1 .* (image1 > 0) );
image2 = image2 - max(image2(:) ) / 300;
image2 = sqrt(image2 .* (image2 > 0) );
image3 = image3 - max(image3(:) ) / 300;
image3 = sqrt(image3 .* (image3 > 0) );

mesh = pixelsize/resize;
[ly, lx] = size(image0);
x = (1:lx)*mesh ; % x = x-mean(x);
y = (1:ly)*mesh;  % y = y-mean(y);
[xx, yy] = meshgrid(x,y);
[~, ii] = max(image0(:) );
x0 = xx(ii);
y0 = yy(ii);
L = (z3 - z0);
xx = gpuArray(xx);
yy = gpuArray(yy);
rxy2 = (xx - x0) .^ 2 + (yy - y0) .^ 2;
iniPhi = k * (sqrt( (rxy2 + L^2) ) - L);  %initial phase
clear image0
clear xx
clear yy
clear rxy2

kmeshx = 2 * pi / (mesh * lx);
kmeshy = 2 * pi / (mesh * ly);
kx = ( (1 : lx) - lx / 2 - 1) * kmeshx;  %kx = kx-mean(kx);
ky = ( (1 : ly) - ly / 2 - 1) * kmeshy;  %ky = ky-mean(ky);
kx = fftshift(kx);
ky = fftshift(ky);
[kkx, kky] = meshgrid(kx, ky);
kkx = gpuArray(kkx);
kky = gpuArray(kky);
kwindow = exp(- (kkx .^ 2 + kky .^ 2) / k^2 / NAs^2);
kwindow = kwindow > (max(kwindow(:) ) / 2.71828);

if(exist([resPath 'currentPhase.mat'], 'file') )
    load([resPath 'currentPhase.mat']);
else
    curPhi = iniPhi;  %current phase
end

gpuImage1 = gpuArray(image1);
gpuImage2 = gpuArray(image2);
gpuImage3 = gpuArray(image3);

if (~(exist('j', 'var') ) )
    j = 0;
end
inij = j;

while j < 100001
    j = j + 1
    % cep stands for current estimated picture
    cep = ifft2(kwindow .* exp(1i * sqrt(k^2 - kkx .^ 2 - kky .^ 2) * (z2 - z3) ) .* fft2(gpuImage3 .* exp(1i * curPhi) ));
    curPhi = angle(cep);
    
    cep = ifft2(kwindow .* exp(1i * sqrt(k^2 - kkx .^ 2 - kky .^ 2) * (z1 - z2) ) .* fft2(gpuImage2 .* exp(1i * curPhi) ));
    curPhi = angle(cep);
    
    cep = ifft2(kwindow .* exp(1i * sqrt(k^2 - kkx .^ 2 - kky .^ 2) * (z3 - z1) ) .* fft2(gpuImage1 .* exp(1i * curPhi) ));
    curPhi = angle(cep);
    
    if (mod(j - inij, 100) == 1)
        subplot(2, 2, 1);
        imagesc(x, y, gpuImage3);  colorbar; title(['iteration = ' num2str(j) ',image3 to image 3']);
        %     subplot(1,3,2);
        %     imagesc(res32); caxis([0 1]);  colorbar; title(['iteration = ' num2str(j) ',image2 to image 3']);
        subplot(2, 2, 2);
        imagesc(x, y, (abs(cep) ) ); colorbar; title(['iteration = ' num2str(j) ',ep3 ']);
        
        subplot(2, 2, 3);
        imagesc(x, y, sin(curPhi - iniPhi) ); colorbar; title(['iteration = ' num2str(j) ',angle ']);
        
        res13 = abs( (abs(cep) ) .^ 2 - (gpuImage3) .^ 2) ./ max( (gpuImage3(:) ) .^ 2);
        res1m( ceil( (j - inij) / 100 ) ) = mean(mean(res13) );
        
        subplot(2,2,4);
        plot(inij+1: 100 : j ,res1m(1, ceil( (j - inij) / 100 ) ),'-');
        drawnow;
        save([resPath 'currentphase.mat'],'curPhi','j');
        %localcurphi = gather(curphi);
    end
end











