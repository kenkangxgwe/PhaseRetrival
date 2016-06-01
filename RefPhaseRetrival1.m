function RefPhaseRetrival1(image0, image1, image2, image3, z0, z1, z2, z3, params, FilePath)
% image0 focal
% image1 small
% image2 medium
% image3 big, but not being clipped by the camera
% after we get the phase for image3, we propagate E3 to the plane where we
% are going to do the experiment.
% An optimization is applied to readjusted the zi.

gpuDevice(1);
wavelength =  params.wavelength;
pixelsize = params.pixelsize;
NAs = params.NAs;
resize = 6;
z10 = z1; z20 = z2; z30 = z3;
k = 2 * pi / wavelength;  % Wave Vector
k0=k;
displayN = 100; % Display the intermediate result every displayN iterations
optimizeN = 500; % Optimize the distance from zi to zi0 every optimizeN iterations

resPath = [FilePath 'Result\'];

image0 = sqrt(double(imresize(image0, resize, 'nearest') ) );

mesh = pixelsize/resize;
[ly, lx] = size(image0);
x = (1 : lx) * mesh; y = (1 : ly) * mesh;
% x = x-mean(x); y = y-mean(y);
[xx, yy] = meshgrid(x,y);
[~, ii] = max(image0(:) );
x0 = xx(ii); y0 = yy(ii);
L = (z3 - z0);

xx = gpuArray(xx); yy = gpuArray(yy);
rxy2 = (xx - x0) .^ 2 + (yy - y0) .^ 2;
iniPhi = k * (sqrt( (rxy2 + L^2) ) - L);  % initial phase
clear image0 xx yy rxy2;

kmeshx = 2 * pi / (mesh * lx); kmeshy = 2 * pi / (mesh * ly);
kx = ( (1 : lx) - lx / 2 - 1) * kmeshx;
ky = ( (1 : ly) - ly / 2 - 1) * kmeshy;
% kx = kx-mean(kx); ky = ky-mean(ky);

kx = fftshift(kx); ky = fftshift(ky);
[kkx, kky] = meshgrid(kx, ky);
kkx = gpuArray(kkx); kky = gpuArray(kky);
kwindow = exp(- (kkx .^ 2 + kky .^ 2) / k^2 / NAs^2);
% kwindow2 = kwindow > (max(kwindow(:) ) / 1.01);
kwindow = kwindow > (max(kwindow(:) ) / 2.71828);
coef1=[0 1]; coef2=[0 1]; coef3=[0 1]; 

if(exist([resPath 'currentPhase.mat'], 'file') )
    load([resPath 'currentPhase.mat']);
else
    curPhi = iniPhi;  % current phase
end

gpuImage1 = real(sqrt(complex(imresize(gpuArray(image1), resize) ) ) );
gpuImage2 = real(sqrt(complex(imresize(gpuArray(image2), resize) ) ) );
gpuImage3 = real(sqrt(complex(imresize(gpuArray(image3), resize) ) ) );
% gpuImage1 = real(ifft2(kwindow2 .* fft2(gpuImage1) ) );
% gpuImage2 = real(ifft2(kwindow2 .* fft2(gpuImage2) ) );
% gpuImage3 = real(ifft2(kwindow2 .* fft2(gpuImage3) ) ); 

if (~(exist('j', 'var') ) )
    j = 0;
end
inij = j;

% Options for fiminc in the iteration
options = optimoptions(@fminunc, 'Display', 'iter', 'Algorithm', 'quasi-newton');
options.TolX = 1e-18; options.TolFun = 1e-18; options.MaxFunEvals = 1000;

z1 = z10 + coef2(1); z2 = z20 + coef1(1); z3 = z30 + coef3(1);

while j < 240000
    j = j + 1
    %     cep stands for current estimated picture
    %      cep = ifft2(kwindow .* exp(1i * sqrt(k^2 - kkx .^ 2 - kky .^ 2) * (z2 - z3) ) .* fft2(gpuImage3 .* exp(1i * curPhi) ) );
    %      curPhi = angle(cep);
    %
    %     cep = ifft2(kwindow .* exp(1i * sqrt(k^2 - kkx .^ 2 - kky .^ 2) * (z1 - z2) ) .* fft2(gpuImage2 .* exp(1i * curPhi) ) );
    %    curPhi = angle(cep);
    
    cep=fft2(gpuImage3 .* exp(1i * curPhi));
    if (mod(j, optimizeN) == 1)
        fun1 = @(x) gather(sum(sum(abs(x(2) * (abs(ifft2(kwindow .* exp(1i * sqrt(k^2 - kkx .^ 2 - kky .^ 2) * (z20 + x(1) - z3 ) ) .* cep) ) ).^2 - gpuImage2.^2) ) ) );
        [coef1] = fminunc(fun1, coef1, options)
        z2=z20 + coef1(1);
    end
    cep = ifft2(kwindow .* exp(1i * sqrt(k^2 - kkx .^ 2 - kky .^ 2) * (z2 - z3) ) .* cep);
    curPhi = angle(cep);
    
    cep=fft2(gpuImage2 .* exp(1i * curPhi));
    if (mod(j, optimizeN) == 2)
        fun1 = @(x) gather(sum(sum(abs(x(2) * (abs(ifft2(kwindow .* exp(1i * sqrt(k^2 - kkx .^ 2 - kky .^ 2) * (z10 + x(1) - z2 ) ) .* cep) ) ).^2 - gpuImage1.^2) ) ) );
        [coef2] = fminunc(fun1, coef2, options)
        z1 = z10 + coef2(1);
    end
    cep=ifft2(kwindow .* exp(1i * sqrt(k^2 - kkx .^ 2 - kky .^ 2) * (z1 - z2 ) ).*cep);
    curPhi = angle(cep);
    
    cep=fft2(gpuImage1 .* exp(1i * curPhi));
    if (mod(j, optimizeN) ==3)
        fun1 = @(x) gather(sum(sum(abs(x(2) * (abs(ifft2(kwindow .* exp(1i * sqrt(k^2 - kkx .^ 2 - kky .^ 2) * (z30 + x(1) - z1 ) ) .* cep) ) ).^2 - gpuImage3.^2) ) ) );
        [coef3] = fminunc(fun1, coef3, options)
        z3 = z30 + coef3(1);
    end
    cep=ifft2(kwindow .* exp(1i *sqrt(k^2 - kkx .^ 2 - kky .^ 2) * (z3 - z1 ) ).*cep);
    curPhi = angle(cep);
    
    if (mod(j, displayN) == 0)
        subplot(2, 2, 1); imagesc(x, y, gpuImage3); colorbar; title(['iteration = ' num2str(j) ', original image3']);        
        subplot(2, 2, 2); imagesc(x, y, coef3(2)^.5 * (abs(cep) ) ); colorbar; title(['iteration = ' num2str(j) ', ep3']);        
        subplot(2, 2, 3); imagesc(x, y, sin(curPhi - iniPhi) ); colorbar; title(['iteration = ' num2str(j) ', angle']);
        % This is a magic function :)
        % residue = abs( coef3(2) * (abs(cep) ).^2  - (gpuImage3).^2) / max( (gpuImage3(:) .^2) );
        a = coef3(2) * (abs(cep) ).^2;
        b = gpuImage3.^2;
        a=a(1500:5100,2500:6300);
        b=b(1500:5100,2500:6300);
        a = a(:); b = b(:);
        % This one calculates the normalized cross correlation
        % residue = 1 - dot(a - mean(a), b - mean(b) ) / (size([a; b], 1) * std(a) * std(b) );
        % This one calculates the cosine of the general angle between two column vectors
        residue = 1 - sum(a.*b) /sqrt(sum(a.^2)*sum(b.^2));
        resArray(ceil( (j - inij) / displayN ) ) = residue;
        subplot(2, 2, 4); plot(inij + 1: displayN : j, resArray(1 : ceil( (j - inij) / displayN ) ), 'x');
        drawnow;
        save([resPath 'currentphase.mat'], 'curPhi', 'j', 'coef1', 'coef2', 'coef3');
    end
end

% Smoothed Phase
delphi = gpuImage3 .* (cos(curPhi - iniPhi) + 1i * sin(curPhi - iniPhi) );
% Averaging filter
ones(10)
curPhiSmoothed = gather(angle(filter2(fspecial('average',10), delphi) ) ) + iniPhi;
save([resPath 'currentphase.mat'], 'curPhiSmoothed', 'iniPhi', '-append');
