clear all

Data=dlmread('I:\170123_TiSa Spectrum\170123_EM2 light source after MMF_7mW.txt');
PSF=Data(:,2)-mean(Data(:,2));

FFT=fft(PSF);
FFT_f=FFT;
FFT_f(1:10)=0;
FFT_f(40:end)=0;
PSF_f=ifft(FFT_f);

plot(Data(:,1),PSF_f,Data(:,1),abs(PSF_f));
plot(Data(:,1),PSF_f);

%%
Pix=0.74/2/58.7143;

Pos=Data(:,1)*Pix;

Env=abs(PSF_f)/max(abs(PSF_f));

FWHM=abs(Pos(find(Env>0.5,1,'first'))-Pos(find(Env>0.5,1,'last')));
figure('Position', [100, 100, 800, 500]);
plot(Pos,PSF_f/max(abs(PSF_f)),Pos,abs(PSF_f)/max(abs(PSF_f)),'LineWidth',1.5);
xlabel('Axial Position (micron)','fontsize',15);
ylabel('Normalzied Power (abs.)','fontsize',15);
xlim([0 11])
legend(sprintf('FWHM=%.03g\\mum',FWHM));
set(gca,'fontsize',15)
