clear all
%%
NA=0.514;

f_Number=1/(2*NA);

Lambda=0.78;    %micron

Rayleigh_Resolution=0.61*Lambda/NA

(1/Rayleigh_Resolution)*1000/2

f0=1/Lambda/f_Number;   %unit: 1/micron


f0_lpm=f0*1000;   %unit: 1/micron

f=0:(f0_lpm*1.1);

index_f_at_Rayleigh=find(f>(1/Rayleigh_Resolution)*1000/2,1,'first');

MTF_dl=2/pi*(acos(f/f0_lpm)-f./f0_lpm.*(1-(f/f0_lpm).^2).^0.5);

MTF_at_Rayleigh=MTF_dl(index_f_at_Rayleigh)


plot(f,MTF_dl);
title(sprintf('Matlab calculated, NA=%g',NA))
ylabel('Spatial Frequency in cycles/mm');
xlabel('Modulus of the OTF');