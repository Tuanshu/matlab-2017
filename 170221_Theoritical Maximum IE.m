clear all

R1=0.29;%

R2=0.15;%

R_Sample=0.04;

T1=(1-R1)*(1-R2)*R_Sample*(1-R2)*(1-R1);
T2=(1-R1)*R2*R1*R2*(1-R1);

DC=T1^2+T2^2;

Inter=2*(T1*T2);

IE=Inter/DC