clear all

x=1:0.01:1000;
a=1;
d=0;
b=500;
c=50;
Function=a./(exp(((x-b)/c)-1))+d;
plot(Function)
Function2=a./(exp(((x-b)/c))+1)+d;1
plot(x,Function2)

`   Function2_diff=abs([diff(Function2) 0]);
PSF=Function2_diff/max(Function2_diff);
plot(x,PSF)


FWHM=abs(x(find(PSF>0.5,1,'first'))-x(find(PSF>0.5,1,'last')))


FWHM/c