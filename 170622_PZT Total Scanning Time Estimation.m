clear all
%%
Multi=64;

a0=0;

a1=0.004184;

a2=2.8e-05;

Time=0:350000;


V=a0+a1*Time.^0.5+a2*Time;
Total_Time=Time(find(V>10,1,'first'))/1000/Multi
plot(Time/1000/Multi,V);
xlim([0 Total_Time]);
