function RefPhaseRetrival1c(image0, image1, image2, image3, z0, z1, z2, z3, pixelsize, wavelength, FilePath)
% In this version, the image is readjusted at each image plane along x-y
% directions. the results seem unstable...
% image0 focal
% image1 small
% image2 medium
% image3 big, but not being clipped by the camera
% after we get the phase for image3, we propagate E3 to the plane where we
% are going to do the experiment.

gpuDevice(1);
resize = 6;
k = 2 * pi / wavelength;  % Wave Vector
NAs = 0.35;
displayN=100; % Display the intermediate result every N iterations

resPath = [FilePath '\Result\'];

image0 = double(imresize(image0, resize, 'nearest') );
image1 = double(imresize(image1, resize, 'nearest') );
image2 = double(imresize(image2, resize, 'nearest') );
image3 = double(imresize(image3, resize, 'nearest') );

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
iniPhi = k * (sqrt( (rxy2 + L^2) ) - L);  % initial phase
clear image0
clear xx
clear yy
clear rxy2

kmeshx = 2 * pi / (mesh * lx);
kmeshy = 2 * pi / (mesh * ly);
kx = ( (1 : lx) - lx / 2 - 1) * kmeshx;  % kx = kx-mean(kx);
ky = ( (1 : ly) - ly / 2 - 1) * kmeshy;  % ky = ky-mean(ky);

kx = fftshift(kx);
ky = fftshift(ky);
[kkx, kky] = meshgrid(kx, ky);
kkx = gpuArray(kkx);
kky = gpuArray(kky);
kwindow = exp(- (kkx .^ 2 + kky .^ 2) / k^2 / NAs^2);
kwindow2 = kwindow > (max(kwindow(:) ) / 1.01);
kwindow = kwindow > (max(kwindow(:) ) / 2.71828);

coef=[0 0 0 1];

if(exist([resPath 'currentPhase.mat'], 'file') )
    load([resPath 'currentPhase.mat']);
else
    curPhi = iniPhi;  % current phase
end

gpuImage1=real(ifft2(kwindow2 .* fft2(gpuArray(image1))));
gpuImage2=real(ifft2(kwindow2 .* fft2(gpuArray(image2))));
gpuImage3=real(ifft2(kwindow2 .* fft2(gpuArray(image3))));
clear kwindow2;

if (~(exist('j', 'var') ) )
    j = 0;
end
inij = j;

% Options for fiminc in the iteration
options = optimoptions(@fminunc,'Display','iter','Algorithm','quasi-newton');
options.TolX=1e-8;
options.TolFun=1e-8;
options.MaxFunEvals= 1000;

while j < 100001
    j = j + 1
%     cep stands for current estimated picture
%     cep = ifft2(kwindow .* exp(1i * sqrt(k^2 - kkx .^ 2 - kky .^ 2) * (z2 - z3) ) .* fft2(gpuImage3 .* exp(1i * curPhi) ) );
%     curPhi = angle(cep);
    
    cep = ifft2(kwindow .* exp(1i * sqrt(k^2 - kkx .^ 2 - kky .^ 2) * (z1 - z3) ) .* fft2(gpuImage3 .* exp(1i * curPhi) ) );
    curPhi = angle(cep);
    
    if (mod(j , displayN) ==0)
        cep=fft2(gpuImage1 .* exp(1i * curPhi));
        fun1 = @(x)gather(sum(sum( abs( (abs(ifft2(kwindow .* exp(1i * sqrt(k^2 - kkx .^ 2 - kky .^ 2) * (z3 - z1+x(1) )+x(2)*kkx.^2+x(3)*kky.^2 ).*cep) ) ).^2*x(4)-gpuImage3.^2) ) ) );
        [coef] = fminunc(fun1, coef, options)
        cep=ifft2(kwindow .* exp(1i * sqrt(k^2 - kkx .^ 2 - kky .^ 2) * (z3 - z1+coef(1) ) +coef(2)*kkx.^2+coef(3)*kky.^2).*cep);
    else
        cep = ifft2(kwindow .* exp(1i * sqrt(k^2 - kkx .^ 2 - kky .^ 2) * (z3 - z1) ) .* fft2(gpuImage1 .* exp(1i * curPhi) ) );
    end
    curPhi = angle(cep);
    
    if (mod(j , displayN) == 0)
        subplot(2, 2, 1);
        imagesc(x, y, gpuImage3);  colorbar; title(['iteration = ' num2str(j) ',image3 to image 3']);
        
        subplot(2, 2, 2);
        imagesc(x, y, coef(4)^.5 * (abs(cep) ) ); colorbar; title(['iteration = ' num2str(j) ',ep3 ']);
        
        subplot(2, 2, 3);
        imagesc(x, y, sin(curPhi - iniPhi) ); colorbar; title(['iteration = ' num2str(j) ',angle ']);
        
        res13 = abs( coef(4) * (abs(cep) ).^2  - (gpuImage3).^2) / max( (gpuImage3(:) .^2) );
        res1m( ceil( (j - inij) / displayN ) ) = mean(mean(res13) );
        
        subplot(2,2,4);
        plot(inij+1: displayN : j ,res1m(1:ceil( (j - inij) / displayN ) ),'x');
        drawnow;
        save([resPath 'currentphase.mat'],'curPhi','j','coef');
    end
end











