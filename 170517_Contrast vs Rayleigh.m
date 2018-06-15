clear all
%%
NA=0.35;

Testing_Line_per_mm=500;

Line_period_micron=1/Testing_Line_per_mm*1000;
Number_of_period=15.625;

Window_Size_Multiplier=1;

f_Number=1/(2*NA);

Lambda=0.55;    %micron

Rayleigh_Resolution=0.61*Lambda/NA;

Sampling_Resolution=0.01;  %micron

X=0:Sampling_Resolution:(Number_of_period*Line_period_micron);

r=(-1*Window_Size_Multiplier*Rayleigh_Resolution):Sampling_Resolution:(Window_Size_Multiplier*Rayleigh_Resolution);
Airy_Profile=real(2*besselj(1,2*pi*(r/(Lambda/NA)))./ (2*pi*(r/(Lambda/NA)))).^2;
%Airy_Profile=gaussmf(r,[0.42*Lambda*f_Number 0]);   %Gauss approx
Airy_Profile(isnan(Airy_Profile))=1;
plot(r,Airy_Profile);
%%

Airy_Profile_Smooth=smooth(Airy_Profile,length(Airy_Profile)/2);

plot(Airy_Profile_Smooth);

%%
Airy_Profile_Shift=circshift(Airy_Profile,round(length(Airy_Profile)/2),2);

plot(r,Airy_Profile,r,Airy_Profile_Shift);

Manual_Sum_Profile=Airy_Profile+Airy_Profile_Shift;

plot(r,Manual_Sum_Profile);

Manual_Contrast=(max(Manual_Sum_Profile)-min(Manual_Sum_Profile))/(max(Manual_Sum_Profile)+min(Manual_Sum_Profile));

Airy_Profile_norm=Airy_Profile/sum(Airy_Profile);

Line_Pattern=(square(X*(2*pi)/Line_period_micron)+1)/2;

plot(X,Line_Pattern);

Convolution_Result=conv(Line_Pattern,Airy_Profile_norm,'valid');

X_valid=(0:1:(length(Convolution_Result)-1))*Sampling_Resolution;

plot(X_valid,Convolution_Result);

Contrast=(max(Convolution_Result)-min(Convolution_Result))/(max(Convolution_Result)+min(Convolution_Result));

