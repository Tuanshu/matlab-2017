clf
clear all
%%
NA=0.51;

f_Number=1/(2*NA);

Lambda=0.78;    %micron

f0=1/Lambda/f_Number;   %unit: 1/micron

Sampling_Resolution=1.35;

Sampling_Freq=1/Sampling_Resolution*1000;

Nyquist=Sampling_Freq/2;

f0_lpm=f0*1000;   %unit: 1/mm

f=0:2000;

MTF_dl=2/pi*(acos(f/f0_lpm)-f./f0_lpm.*(1-(f/f0_lpm).^2).^0.5);

Rayleigh=0.61*Lambda/NA;

hax=axes; 
plot(f,MTF_dl);
title('Lateral Resolution')
hold on
line([Nyquist Nyquist],get(hax,'XLim'),'LineWidth',2,'Color','green')
xlim([0 Sampling_Freq*2]);
ylim([0 1]);
xlabel('Spatial Frequency (line/mm)');
ylabel('SFR (MTF)');
hold off
