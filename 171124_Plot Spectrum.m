Spectrum=dlmread('D:\50perc.txt');

plot(Spectrum(:,1),Spectrum(:,2)/max(Spectrum(:,2)),'LineWidth',2)
xlabel('Wavelength (nm)')
ylabel('Normalized Spectral Density (a.u.)')