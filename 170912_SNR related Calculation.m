clear all
%%
FWC_e=200000;
Bit=12;
DN_max=2^Bit;
Gain=DN_max/FWC_e;
X_number=1024;
Y_number=11;
Sturation_Level=0.9;
PreAVE=2;
PostAVE=2;
Shot_e=((PreAVE*PostAVE*FWC_e*Sturation_Level)^0.5)/PreAVE/PostAVE;
Shot_dn=Shot_e*Gain;