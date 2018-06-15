clear all
fclose all
%%
Camera_Frame_Rate=2501; %fps
Assumed_Central_Wavelength=0.78;    %micron
RI=1.406;
Multi=64;

%% About EMScanner Setting
P_start=0;  %wrong calibration
P_end=290;  %wrong calibration
Total_Number_of_Frame_B_scan_real=5672*2;    %*2 for 2 frames per files
Total_Scanning_Time_real=Total_Number_of_Frame_B_scan_real/Camera_Frame_Rate

V_start=0;  %volt
V_end=10;   %volt
%% About currently input waveform
Time_1stguess=0:1:200000;    %ms

a0_Old=0;
a1_Old=0.004184;
a2_Old=2.8e-05; %actually 2.783E-5

Current_Voltage_Waveform_1stguess=a0_Old+(a1_Old*(Time_1stguess*Multi).^0.5)+a2_Old*Time_1stguess*Multi;
Total_Time=Time_1stguess(find(Current_Voltage_Waveform_1stguess<V_end,1,'last'));         %Always start at 0
Time=Time_1stguess(1:find(Current_Voltage_Waveform_1stguess<V_end,1,'last'));
Current_Voltage_Waveform=Current_Voltage_Waveform_1stguess(1:find(Current_Voltage_Waveform_1stguess<V_end,1,'last'));
plot(Time,Current_Voltage_Waveform);

%% 1st order derivative of V (dVdt)
dVdt_Waveform=a1_Old*Time.^(-0.5)+a2_Old;

plot(Time,dVdt_Waveform);


%% Read BinFPC
Folder_Path='Z:\研究開發部\何端書\AMO\170713_PZT calibration\Input 7\Processed Data\txt\';

file_list=dir(Folder_Path);
Original_Lnegth_file_list=length(file_list);
for p=1:Original_Lnegth_file_list
    if file_list(Original_Lnegth_file_list+1-p).isdir ~= 0
        file_list(Original_Lnegth_file_list+1-p)=[];
    elseif strcmp(file_list(Original_Lnegth_file_list+1-p).name,'.') == 1
        file_list(Original_Lnegth_file_list+1-p)=[];
    elseif strcmp(file_list(Original_Lnegth_file_list+1-p).name,'..') == 1
        file_list(Original_Lnegth_file_list+1-p)=[];
    elseif strcmp(file_list(Original_Lnegth_file_list+1-p).name,'Processed Data') == 1
        file_list(Original_Lnegth_file_list+1-p)=[];
    end
end
BinFPC=0;
for p=1:length(file_list)
    File_Path=[Folder_Path file_list(p).name];
    BinFPC=BinFPC+dlmread(File_Path);
end
BinFPC=BinFPC/length(file_list);
Bin=BinFPC(:,1);
FPC=BinFPC(:,2);

plot(Bin,FPC);
%% Excluding Some Data Point
Bin_Not_Using=[1 2 20 21 22];
FPC(Bin_Not_Using)=NaN;

%% Generating High-resolution FPC waveform
Number_of_Frame_per_Bin=500;    %Time_Window

Total_Time_Recorded=length(Bin)*Number_of_Frame_per_Bin/Camera_Frame_Rate;

Time_Bin=Number_of_Frame_per_Bin/Camera_Frame_Rate*(Bin-0.5)*1000;   %ms, move the time tag to the central instead of the end of bin

FPC_hires=interp1(Time_Bin,FPC,Time);

Time_noNaN=Time;
Time_noNaN(isnan(FPC_hires))=[];
FPC_hires_noNaN=FPC_hires;
FPC_hires_noNaN(isnan(FPC_hires))=[];

plot(Time_Bin,FPC,Time_noNaN,FPC_hires_noNaN);

MSE=sum((8-FPC_hires_noNaN).^2)
%% Position Waveform Generation (Method 1)
Velocity_Waveform_real=Camera_Frame_Rate./FPC_hires_noNaN*Assumed_Central_Wavelength/2/RI;
plot(Time_noNaN,Velocity_Waveform_real);
Position_Waveform=cumsum(Velocity_Waveform_real)/1000;
plot(Time_noNaN,Position_Waveform);

%% Find the corresponding voltage array (for both Method 1&2)
Voltage_Waveform_of_Interest=interp1(Time,Current_Voltage_Waveform,Time_noNaN);
plot(Time_noNaN,Velocity_Waveform_real);


%% Calculate dVdt (Method 2)
dVdt_Waveform_of_Interest=interp1(Time,dVdt_Waveform,Time_noNaN);
plot(dVdt_Waveform_of_Interest,Velocity_Waveform_real);

%% 因起始位置不明, 此方式不能fit a0, 合理 (Method 1)
Expect_Speed=Camera_Frame_Rate/4/2*Assumed_Central_Wavelength/2/RI;
Time_New=Position_Waveform/Expect_Speed*1000+Time_noNaN(1); %Add a Time offset based on the original data
Voltage_Waveform_New=Voltage_Waveform_of_Interest;
plot(Time_New,Voltage_Waveform_New);

%% Generating Equal spacing Array (Method 1)
Voltage_Waveform_New_Equalspacing=interp1(Time_New,Voltage_Waveform_New,Time);

plot(Time_New,Voltage_Waveform_New,Time,Voltage_Waveform_New_Equalspacing);

Voltage_Waveform_New_Equalspacing_noNaN=Voltage_Waveform_New_Equalspacing;
Time_New_NoNaN=Time;
Voltage_Waveform_New_Equalspacing_noNaN(isnan(Voltage_Waveform_New_Equalspacing))=[];
Time_New_NoNaN(isnan(Voltage_Waveform_New_Equalspacing))=[];
Time_New_NoNaN=Time_New_NoNaN*Multi;
plot(Time_New_NoNaN,Voltage_Waveform_New_Equalspacing_noNaN);

%% Fitting (Method 1)

gs = fittype( @(a0,a1, a2, x) a0+a1*x.^0.5+a2*x);
COEFS=fit(Time_New_NoNaN',Voltage_Waveform_New_Equalspacing_noNaN',gs,'StartPoint',[0 0.004 2.8E-5]);%,'StartPoint',[0 0 COEF2.a COEF2.b COEF2.c]);


plot(Time_New_NoNaN,Voltage_Waveform_New_Equalspacing_noNaN,Time_New_NoNaN,COEFS.a0+COEFS.a1*Time_New_NoNaN.^0.5+COEFS.a2*Time_New_NoNaN);
plot(Time,COEFS.a1*(Time*Multi).^0.5+COEFS.a2*(Time*Multi),Time,a0_Old+a1_Old*(Time*Multi).^0.5+a2_Old*(Time*Multi));
ylim([0 10])
%%
a0_New=0

a1_New=COEFS.a1

a2_New=COEFS.a2
