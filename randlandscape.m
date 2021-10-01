function matt = randlandscape(szs,prms,randkey)
% Create a random surfrace using FFT
% matt = randlandscape(szs,prms,randkey)
% szs gives (x,y) size of surface
% prms=[rho,kc] give color and cutoff

if(nargin<3) randkey=0; end;
if(length(szs)<2)
    szs(2)=szs(1);
end;
baseval = 5; % value to be put at (kx,ky)=(0,0)

% two main parameters, choosing color and cutoff
rho=prms(1);
kw =prms(2)*max(szs)/2;

% normalize in case of size-x and size-y not the same
norms=szs/max(szs);

% create a grid to calculate the distance from the center
[xx,yy]=meshgrid(-szs(1)/2:szs(1)/2-1,-szs(2)/2:szs(2)/2-1);
radius = sqrt((xx/norms(1)).^2+(yy/norms(2)).^2);

% power-law with a cutoff enacted on the distance matrix
F=(radius.^-rho).*exp(-radius/kw);

% randomize if needed
if(randkey) 
    rng(randkey);
end;
% create the noise (uniform noise between -1 and 1, per pixel)
noise = rand(size(F))*2-1;
% shift to prepare for the FFT
tmp=fftshift(F.*noise);
tmp(1,1)=baseval;

% Calculate FFT, add real and imaginary values so there's no symmetry
matt=real(fft2(tmp))+imag(fft2(tmp));

end