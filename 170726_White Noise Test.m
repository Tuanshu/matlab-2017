clear all

Noise=wgn(1024,1,0);

Noise_d=diff(Noise);

FFT_Noise=fft(Noise);
FFT_Noise_d=fft(Noise_d);

plot(abs(FFT_Noise));
plot(abs(FFT_Noise_d));
plot(Noise_d);

var(Noise_d)