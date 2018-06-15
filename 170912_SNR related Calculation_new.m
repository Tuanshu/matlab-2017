clear all
%% System Parameter
FWC_e=200000;
Bit=12;
DN_max=2^Bit;
Dark_Noise_dn=0.5*2^4; %DN @12bit
Gain=DN_max/FWC_e;
X_number=1024;
Y_number=11;
Sturation_Level=0.92;
DN_mean=DN_max*Sturation_Level;
PreAVE=2;
Npoint_Corr=0.5;   %experiemental unbiased estimator for N=4
Shot_e=((PreAVE*FWC_e*Sturation_Level)^0.5)/PreAVE;
Shot_dn=Shot_e*Gain;
Shot_PostN_dn=Shot_dn*Npoint_Corr;
IE=0.02;
OPD=5;  %micron
Fringe_Period_frequency=3E8/((OPD*1E-6)*2);

%% Simulation Setting
Repeat_Times=2000;  %for Noise Calculation
%% Spectrum Parameter
X_spectrum_number=1024;
Center_Wavelengh=0.78;  %micron
Wavelength_FWHM=0.17;   %micron
Spectrum_Sampling_Ratio=1;    %?Spectrum?FWHM????????
ZP_Ratio=8;
%% SD-OCT Part
%% Simulating Spectrum (Linear with k, simplified) without Noise
Frequency_FWHM=3E8*Wavelength_FWHM/(Center_Wavelengh^2)/(1E-6);
Calculated_Axial_Resolution=3E8/Frequency_FWHM/2*1E6;   %micron
Center_Frequency=3E8/(Center_Wavelengh*1E-6);
Frequency_Sampling_Resolution=Frequency_FWHM/Spectrum_Sampling_Ratio/X_spectrum_number;

Frequency_Array=Frequency_Sampling_Resolution*((-1*X_spectrum_number/2+1):X_spectrum_number/2)+round(Center_Frequency/Frequency_Sampling_Resolution)*Frequency_Sampling_Resolution;
% DC
Spectrum_e=gaussmf(Frequency_Array,[Frequency_FWHM/2.355 Center_Frequency]).*FWC_e.*Sturation_Level;
Spectrum_dn=Spectrum_e*Gain;    %?????????Ave spectrum?, ????round
% Inter
Inter_Spectrum_e=Spectrum_e*IE.*cos((Frequency_Array-Center_Frequency)/Fringe_Period_frequency*2*pi);
Spectrum_with_Inter_e=Spectrum_e+Inter_Spectrum_e;
% Full Frequency Range
Full_Frequency_Array=[0:Frequency_Sampling_Resolution:ZP_Ratio*max(Frequency_Array)]';
Min_Frequency_Index=find(Full_Frequency_Array>min(Frequency_Array),1,'first');
%% Generating Position Array
Time_max=1/Frequency_Sampling_Resolution;    %sec
Temporal_Resolution=Time_max/length(Full_Frequency_Array);
Time=[Temporal_Resolution:Temporal_Resolution:Time_max];
Position_Full=3E8*Time/2;  %/2 for round-trip
Position=Position_Full(1:round(size(Full_Frequency_Array)/2));
Position_micron=Position*1E6;
Position_Sampling_Resolution_micron=Position_micron(2)-Position_micron(1);
Index_Theoritical_Peak_Signal=round(OPD/Position_Sampling_Resolution_micron);

%% Simulating Spectrum (Linear with k, simplified) with Noise
Signal_SDOCT_Array=zeros(Repeat_Times,1);
for p=1:Repeat_Times
    % Noise
    Shot_Noise_Spectrum_e=randn(X_spectrum_number,1)'.*Spectrum_with_Inter_e.^0.5;

    Spectrum_e_with_Noise=Spectrum_with_Inter_e+Shot_Noise_Spectrum_e;
    Spectrum_dn_with_Noise=round(Spectrum_e_with_Noise.*Gain);
    % plot(Frequency_Array,Spectrum_dn_with_Noise);
    %% SD-OCT Calculation (Tranditional)
    Full_Spectrum=zeros(size(Full_Frequency_Array,1),size(Full_Frequency_Array,2));
    Full_Spectrum(Min_Frequency_Index:(Min_Frequency_Index+length(Spectrum_dn_with_Noise)-1))=Spectrum_dn_with_Noise-Spectrum_dn;
    % plot(Full_Frequency_Array,Full_Spectrum);

    S_Full=fft(Full_Spectrum);
    S=S_Full(1:round(size(S_Full)/2));
    % 
    % plot(abs(S))
    % xlim([50 300]/8*ZP_Ratio)
    % ylim([-4E4 4E4])
    
    % Record Signal
    Signal_SDOCT_Array(p)=abs(S(Index_Theoritical_Peak_Signal));
    disp(p);
end
%% Plot Signal

plot(Position_micron,abs(S),Position_micron,real(S))
xlim([OPD-5 OPD+10])
ylim([-4E4 4E4])

%% SNR Calculation for SD-OCT
plot(Signal_SDOCT_Array);

Signal_SDOCT=mean(Signal_SDOCT_Array);
Noise_SDOCT=std(Signal_SDOCT_Array);
SNR_SDOCT=20*log10(Signal_SDOCT/Noise_SDOCT);

%% TD-OCT Part (No averaging and No N-point)
%% Simulating Axial Signal Profile for TD-OCT without Noise
%% (Note about the bandpass issue)
Frame_per_Carrier=8;
Axial_Length=400;        %micron
Position_of_Inter=5;    %micron
% Inter (1st)
%% Generating TD Axial Profile (All unit of length is regarding to OPD, not actual thickness)
Axial_Sampling_Resolution_micron=Center_Wavelengh/Frame_per_Carrier;
Axial_Position_micron=Axial_Sampling_Resolution_micron:Axial_Sampling_Resolution_micron:round(Axial_Length/Axial_Sampling_Resolution_micron)*Axial_Sampling_Resolution_micron;
Axial_DC_e=ones(1,length(Axial_Position_micron)).*FWC_e.*Sturation_Level;
Axial_Envelope_e=IE*gaussmf(Axial_Position_micron,[Calculated_Axial_Resolution/2.355 Position_of_Inter]).*FWC_e.*Sturation_Level;
Axial_Carrier_norm=cos((Axial_Position_micron-Position_of_Inter)/(Center_Wavelengh/2)*2*pi);
Axial_Inter_e=Axial_Envelope_e.*Axial_Carrier_norm;
Axial_Signal_e=Axial_DC_e+Axial_Inter_e;
Axial_Inter_dn=Axial_Inter_e*Gain;  %Note, ??round, ????????SNR?
[value Index_Peak_Signal]=min(abs(Axial_Position_micron-Position_of_Inter));

%% TD-OCT SNR Calculation
Signal_TDOCT_Array=zeros(Repeat_Times,1);
for p=1:Repeat_Times
    Axial_Shot_Noise_e=randn(length(Axial_Signal_e),1)'.*(Axial_Signal_e.^0.5);
    Axial_Signal_with_noise_e=Axial_Signal_e+Axial_Shot_Noise_e;
    Axial_Signal_with_noise_dn=round(Axial_Signal_with_noise_e*Gain);
    %plot(Axial_Position_micron,Axial_Signal_with_noise_dn);
    Signal_TDOCT_Array(p)=Axial_Signal_with_noise_dn(Index_Peak_Signal);
    disp(p)
end

Signal_TDOCT=Axial_Inter_dn(Index_Peak_Signal);       %mean(Signal_TDOCT_Array)?????DC?
Noise_TDOCT=std(Signal_TDOCT_Array);
SNR_TDOCT=20*log10(Signal_TDOCT/Noise_TDOCT);

%% SNR between TD-OCT and SD-OCT

SNR_diff=SNR_SDOCT-SNR_TDOCT

Theoritical_SNR_diff=10*log10(X_spectrum_number/2)

Ratio=10^(SNR_diff/20);

%% Speed Related Calculation
% Calculate the number of electrons required for above TD and SD schme
Total_Electron_SDOCT=sum(Spectrum_e)*X_number;
Total_Electron_TDOCT=FWC_e.*Sturation_Level*X_number;

FrameRate_TDOCT=4000;   %Partial Frame, 1024*11, fps
FrameRate_SDOCT=FrameRate_TDOCT/Total_Electron_SDOCT*Total_Electron_TDOCT;    %Full Frame, 1024*1024, fps

Depth_per_Frame_SDOCT=Center_Wavelengh/2;    %micron
Depth_per_Frame_TDOCT=Center_Wavelengh/2;   %micron

Target_Depth=400;   %micron

Frame_Required_for_Target_Depth_SDOCT=Target_Depth/Depth_per_Frame_SDOCT;
Frame_Required_for_Target_Depth_TDOCT=Target_Depth/Depth_per_Frame_TDOCT;

Time_Required_for_Target_Depth_SDOCT=Frame_Required_for_Target_Depth_SDOCT/FrameRate_SDOCT;
Time_Required_for_Target_Depth_TDOCT=Frame_Required_for_Target_Depth_TDOCT/FrameRate_TDOCT;

Avaliable_Averaging_TDOCT=round(Time_Required_for_Target_Depth_SDOCT/Time_Required_for_Target_Depth_TDOCT);
% Ave Factor = 8; ?????Include 4-point and PreAVE???(????PostAVE, ???Depth
% per Frame??????PostAVE)

%%

%% TD-OCT Part (with AVE)
%% Generating TD Axial Profile (All unit of length is regarding to OPD, not actual thickness)
Axial_DC_AVE_e=ones(1,length(Axial_Position_micron)).*FWC_e.*Sturation_Level*Avaliable_Averaging_TDOCT;
Axial_DC_AVE_dn=Axial_DC_AVE_e*Gain;    %?round, ????????
Axial_Envelope_AVE_e=IE*gaussmf(Axial_Position_micron,[Calculated_Axial_Resolution/2.355 Position_of_Inter]).*FWC_e.*Sturation_Level*Avaliable_Averaging_TDOCT;
Axial_Inter_AVE_e=Axial_Envelope_AVE_e.*Axial_Carrier_norm;
Axial_Signal_AVE_e=Axial_DC_AVE_e+Axial_Inter_AVE_e;
Axial_Inter_AVE_dn=Axial_Inter_AVE_e*Gain;  %Note, ??round, ????????SNR?

%% TD-OCT (AVE) SNR Calculation
Signal_TDOCT_AVE_Array=zeros(Repeat_Times,1);
for p=1:Repeat_Times
    Axial_Shot_Noise_AVE_e=randn(length(Axial_Signal_AVE_e),1)'.*(Axial_Signal_AVE_e.^0.5);
    Axial_Signal_with_noise_AVE_e=Axial_Signal_AVE_e+Axial_Shot_Noise_AVE_e;
    Axial_Signal_with_noise_AVE_dn=Axial_Signal_with_noise_AVE_e*Gain; %???AVE, ???round?
    %plot(Axial_Position_micron,Axial_Signal_with_noise_dn);
    Signal_TDOCT_AVE_Array(p)=Axial_Signal_with_noise_AVE_dn(Index_Peak_Signal);
    disp(p)
end

Signal_TDOCT_AVE=Axial_Inter_AVE_dn(Index_Peak_Signal);       %mean(Signal_TDOCT_Array)?????DC?
Noise_TDOCT_AVE=std(Signal_TDOCT_AVE_Array);
SNR_TDOCT_AVE=20*log10(Signal_TDOCT_AVE/Noise_TDOCT_AVE);


%% TD-OCT Part (with AVE and Bandpass Filtering)
%% TD-OCT (AVE + BF) SNR Calculation
Signal_TDOCT_AVE_BF_Array=zeros(Repeat_Times,1);
for p=1:Repeat_Times
    Axial_Shot_Noise_AVE_BF_e=randn(length(Axial_Signal_AVE_e),1)'.*(Axial_Signal_AVE_e.^0.5);
    Axial_Signal_with_noise_AVE_BF_e=Axial_Signal_AVE_e+Axial_Shot_Noise_AVE_BF_e;
    Axial_Signal_with_noise_AVE_BF_dn=Axial_Signal_with_noise_AVE_BF_e*Gain; %???AVE, ???round?
    %plot(Axial_Position_micron,Axial_Signal_with_noise_dn);
    %% Bandpass Filtering
    %plot(real(fft(Axial_Inter_AVE_dn)))
    Axial_Signal_with_noise_AVE_BF_DCremoved_dn=Axial_Signal_with_noise_AVE_BF_dn-Axial_DC_AVE_dn;
    FFT=fft(Axial_Signal_with_noise_AVE_BF_DCremoved_dn);
    RP=725;
    LP=1315;
    RP_c=length(FFT)-RP;
    LP_c=length(FFT)-LP;
    
    FFT_filtered=FFT;
    FFT_filtered(1:RP)=0;
    FFT_filtered(LP:LP_c)=0;
    FFT_filtered(RP_c:end)=0;
    Axial_Signal_with_noise_AVE_BF_Filtered_dn=real(ifft(FFT_filtered));
    %plot(Axial_Position_micron,Axial_Signal_with_noise_AVE_BF_Filtered_dn,Axial_Position_micron,Axial_Signal_with_noise_AVE_BF_DCremoved_dn)
    Signal_TDOCT_AVE_BF_Array(p)=Axial_Signal_with_noise_AVE_BF_Filtered_dn(Index_Peak_Signal);
    disp(p)
end

Signal_TDOCT_AVE_BF=mean(Signal_TDOCT_AVE_BF_Array);    %Signal only, ??????mean, ???????
Noise_TDOCT_AVE_BF=std(Signal_TDOCT_AVE_BF_Array);
SNR_TDOCT_AVE_BF=20*log10(Signal_TDOCT_AVE_BF/Noise_TDOCT_AVE_BF);


%% SD-OCT (Roll-off optimization, reducing Spectral Resolution)
Spectral_AVE_Factor=4;
Frequency_Sampling_Resolution_SRreduced=Frequency_Sampling_Resolution*Spectral_AVE_Factor;
Full_Frequency_SRreduced_Array=[0:Frequency_Sampling_Resolution_SRreduced:ZP_Ratio*max(Frequency_Array)]';
Min_Frequency_SRreduced_Index=find(Full_Frequency_SRreduced_Array>min(Frequency_Array),1,'first');

%% Generating Position Array (Roll-offed)
Time_max_SRreduced=1/Frequency_Sampling_Resolution_SRreduced;    %sec
Temporal_Resolution_SRreduced=Time_max_SRreduced/length(Full_Frequency_SRreduced_Array);%same
Time_SRreduced=[Temporal_Resolution_SRreduced:Temporal_Resolution_SRreduced:Time_max_SRreduced];
Position_Full_SRreduced=3E8*Time_SRreduced/2;  %/2 for round-trip
Position_SRreduced=Position_Full_SRreduced(1:round(size(Full_Frequency_SRreduced_Array)/2));
Position_micron_SRreduced=Position_SRreduced*1E6;
Position_Sampling_Resolution_micron_SRreduced=Position_micron_SRreduced(2)-Position_micron_SRreduced(1);
Index_Theoritical_Peak_Signal_SRreduced=round(OPD/Position_Sampling_Resolution_micron_SRreduced);

Signal_SDOCT_SRreduced_Array=zeros(Repeat_Times,1);
for p=1:Repeat_Times
    % Add Noise
    Shot_Noise_Spectrum_e=randn(X_spectrum_number,1)'.*Spectrum_with_Inter_e.^0.5;
    Spectrum_e_with_Noise=Spectrum_with_Inter_e+Shot_Noise_Spectrum_e;
    Spectrum_dn_with_Noise=round(Spectrum_e_with_Noise.*Gain);
    plot(Frequency_Array,Spectrum_dn_with_Noise);
    % Reducing Spectral Resolution
    Reduced_Spectrum_Width=floor(X_spectrum_number/Spectral_AVE_Factor);
    Temp_Spectrum=0;
    Temp_Spectrum_DC=0;
    for q=1:Spectral_AVE_Factor
        Temp_Spectrum=Temp_Spectrum+Spectrum_dn_with_Noise(q:Spectral_AVE_Factor:(Reduced_Spectrum_Width*Spectral_AVE_Factor-Spectral_AVE_Factor+q));
        Temp_Spectrum_DC=Temp_Spectrum_DC+Spectrum_dn(q:Spectral_AVE_Factor:(Reduced_Spectrum_Width*Spectral_AVE_Factor-Spectral_AVE_Factor+q));
    
    end
    Spectrum_dn_with_Noise_SRreduced=Temp_Spectrum/Spectral_AVE_Factor;
    Spectrum_dn_SRreduced=Temp_Spectrum_DC/Spectral_AVE_Factor;

    %% SD-OCT Calculation (Tranditional)
    Full_Spectrum_SRreduced=zeros(size(Full_Frequency_SRreduced_Array,1),size(Full_Frequency_SRreduced_Array,2));
    Full_Spectrum_SRreduced(Min_Frequency_SRreduced_Index:(Min_Frequency_SRreduced_Index+length(Spectrum_dn_with_Noise_SRreduced)-1))=Spectrum_dn_with_Noise_SRreduced-Spectrum_dn_SRreduced;
    % plot(Full_Frequency_Array,Full_Spectrum);

    S_Full_SRreduced=fft(Full_Spectrum_SRreduced);
    S_SRreduced=S_Full_SRreduced(1:round(size(S_Full_SRreduced)/2));
    % 
    % plot(abs(S))
    % xlim([50 300]/8*ZP_Ratio)
    % ylim([-4E4 4E4])
    
    % Record Signal
    Signal_SDOCT_SRreduced_Array(p)=abs(S_SRreduced(Index_Theoritical_Peak_Signal_SRreduced));
    disp(p);
end
%% Plot Signal

plot(Position_micron,abs(S),Position_micron,real(S),Position_micron_SRreduced,abs(S_SRreduced),Position_micron_SRreduced,real(S_SRreduced))
xlim([OPD-5 OPD+10])
ylim([-4E4 4E4])

%% SNR Calculation for SD-OCT (SRreduced)
plot(Signal_SDOCT_SRreduced_Array);

Signal_SDOCT_SRreduced=mean(Signal_SDOCT_SRreduced_Array);
Noise_SDOCT_SRreduced=std(Signal_SDOCT_SRreduced_Array);
SNR_SDOCT_SRreduced=20*log10(Signal_SDOCT_SRreduced/Noise_SDOCT_SRreduced);


%% SD-OCT (Roll-off optimization, Pixel smoothing)

Signal_SDOCT_Smooth_Array=zeros(Repeat_Times,1);
for p=1:Repeat_Times
    % Add Noise
    Shot_Noise_Spectrum_e=randn(X_spectrum_number,1)'.*Spectrum_with_Inter_e.^0.5;
    Spectrum_e_with_Noise=Spectrum_with_Inter_e+Shot_Noise_Spectrum_e;
    Spectrum_dn_with_Noise=round(Spectrum_e_with_Noise.*Gain);
    plot(Frequency_Array,Spectrum_dn_with_Noise);
    % Smoothing
    Spectrum_dn_with_Noise_Smooth=smooth(Spectrum_dn_with_Noise,Spectral_AVE_Factor);
    % SD-OCT Calculation (Tranditional)
    Full_Spectrum_Smooth=zeros(size(Full_Frequency_Array,1),size(Full_Frequency_Array,2));
    Full_Spectrum_Smooth(Min_Frequency_Index:(Min_Frequency_Index+length(Spectrum_dn_with_Noise)-1))=Spectrum_dn_with_Noise-Spectrum_dn;
    % plot(Full_Frequency_Array,Full_Spectrum);

    S_Full_Smooth=fft(Full_Spectrum_Smooth);
    S_Smooth=S_Full_Smooth(1:round(size(S_Full_Smooth)/2));
    % 
    % plot(abs(S))
    % xlim([50 300]/8*ZP_Ratio)
    % ylim([-4E4 4E4])
    
    % Record Signal
    Signal_SDOCT_Smooth_Array(p)=abs(S_Smooth(Index_Theoritical_Peak_Signal));
    disp(p);
end
%%
plot(Position_micron,abs(S_Smooth),Position_micron,real(S_Smooth))
xlim([OPD-5 OPD+10])
ylim([-4E4 4E4])

%% SNR Calculation for SD-OCT (SRreduced)

Signal_SDOCT_Smooth=mean(Signal_SDOCT_Smooth_Array);
Noise_SDOCT_Smooth=std(Signal_SDOCT_Smooth_Array);
SNR_SDOCT_Smooth=20*log10(Signal_SDOCT_Smooth/Noise_SDOCT_Smooth);

%% Final SNR Comparision
SNR_diff_final=SNR_SDOCT_SRreduced-SNR_TDOCT_AVE_BF
Ratio_final=10^(SNR_diff_final/20)
