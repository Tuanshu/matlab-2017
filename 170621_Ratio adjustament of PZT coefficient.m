clear all

Supposed_N=8;

Measured_N=7.5;


Old_a0=0;

Old_a1=0.004321;

Old_a2=2.9687E-5;


New_a0=Old_a0

New_a1=Old_a1*(Measured_N/Supposed_N)^0.5

New_a2=Old_a2*Measured_N/Supposed_N