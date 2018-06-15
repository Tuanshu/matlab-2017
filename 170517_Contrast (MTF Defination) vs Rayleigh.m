clear all
%%
NA=0.45;

f_Number=1/(2*NA);

Lambda=0.55;    %micron

Rayleigh_Resolution=0.61*Lambda/NA;

r=0:0.001:1;
Airy_Profile=real(2*besselj(1,2*pi*r*Lambda/NA)./ (2*pi*r*Lambda/NA)).^2;

plot(r,Airy_Profile);