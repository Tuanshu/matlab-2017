clear all

Supposed_N=8;

Measured_N=8*4000/2500;


Old_a0=0;

Old_a1=0.0048;

Old_a2=0.00002783;

Old_Scanning_Speed=1.367;

New_a0=Old_a0

New_a1=Old_a1*(Measured_N/Supposed_N/2)^0.5

New_a2=Old_a2*Measured_N/Supposed_N/2

New_Scanning_Speed=Old_Scanning_Speed*Measured_N/Supposed_N/2