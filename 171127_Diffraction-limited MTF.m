clf
clear all
%
NA=0.48;

f_Number=1/(2*NA);

Lambda=0.78;    %micron

f0=1/Lambda/f_Number;   %unit: 1/micron

Sampling_Resolution=0.45;

Sampling_Freq=1/Sampling_Resolution*1000;

Nyquist=Sampling_Freq/2;

f0_lpm=f0*1000;   %unit: 1/mm

f=0:3000;

MTF_dl=2/pi*(acos(f/f0_lpm)-f./f0_lpm.*(1-(f/f0_lpm).^2).^0.5);
MTF50=f(find(MTF_dl<0.5,1,'first'));
MTF30=f(find(MTF_dl<0.3,1,'first'));
MTF35=f(find(MTF_dl<0.35,1,'first'));
MTF09=f(find(MTF_dl<0.09,1,'first'));
MTFat500lpsm=MTF_dl(find(f>500,1,'first'))

Rayleigh=0.61*Lambda/NA;
h=plot(f,MTF_dl,'Linewidth',2);
hold on
%line([Nyquist Nyquist],get(hax,'XLim'),'LineWidth',2,'Color','green')
xlim([0 1000]);
ylim([0 1]);
xlabel('Spatial Frequency (line/mm)','FontSize',12);
ylabel('MTF','FontSize',12);
title(sprintf('MTF for Wavelength=%gnm, NA=%g \n(MTF50=%glpm, MTF30=%glpm, MTF35=%glpm, MTF09=%glpm, Rayleigh=%gmicron, MTF@500lps/mm=%g)',1000*Lambda,NA,MTF50,MTF30,MTF35,MTF09,Rayleigh,MTFat500lpsm))
hold off
saveas(h,sprintf('%gnm_NA%g.png',1000*Lambda,NA));