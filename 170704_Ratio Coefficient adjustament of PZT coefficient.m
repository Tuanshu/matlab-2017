clear all
fclose all

Time=0:1:10000;    %ms

Multi=64;

Old_a0=-0.1726;

Old_a1=0.003324049893174;

Old_a2=3.124999999999999e-05;


a0_corr=-1;
a1_corr=0.3;
a2_corr=-0.05;


New_a0=Old_a0*(1+a0_corr)

New_a1=Old_a1*(1+a1_corr)

New_a2=Old_a2*(1+a2_corr)


V=New_a0+New_a1*(Time*Multi).^0.5+New_a2*(Time*Multi);

plot(Time,V);
ylim([0 10]);
xlim([0 Time(find(V>10,1,'first'))]);
xlabel('Time (ms)');
ylabel('V (Voltage)');

Slope=diff(V)/(Time(2)-Time(1));

Slope_2V=Slope(find(V>2,1,'first'))
Slope_5V=Slope(find(V>5,1,'first'))
Slope_8V=Slope(find(V>8,1,'first'))