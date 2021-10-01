function corfunc=spacorr2d(img)
% get the image in the spectral domain
fimg=fft2(img);
% apply convolution by multiplication and inverse-FFT
tmpcf = (ifft2(fimg.* conj(fimg)));
% normalize by maximum value
tmpcf = tmpcf/max(tmpcf(:));
% average the vertical and horizontal axes
corfunc=(tmpcf(1:length(img)/2,1)+tmpcf(1,1:length(img)/2)')/2;

end