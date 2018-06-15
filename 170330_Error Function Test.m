clear all

x=1:1000;

FWHM=100;

Gaussian=gaussmf(x,[FWHM/((log(2)*2)^0.5)/2 500]);

plot(x,Gaussian)

ESF=cumsum(Gaussian);

plot(x,ESF/max(ESF))
