clear all
fclose all
%%
Number_of_Frame_per_Bin=200;    %Time_Window*Ave_PreSpectrogram
If_DifferentModelFor2Seg=1;
%%
Camera_Frame_Rate=4000; %fps
Assumed_Central_Wavelength=0.78;    %micron
RI=1.406;
Old_Multi=128;
Interp_Target_Multi=64;
%% About EMScanner Setting
P_start=0;  %wrong calibration
P_end=400;  %wrong calibration
v1x=1.0936; %micron/sec

%% About currently input waveform
Time_1stguess=0:1:200000;    %ms

a0_Old=10;
a1_Old=-0.0073;
a2_Old=-1.69E-05; %actually 2.783E-5

%%
Total_Time_Calc=abs(P_end-P_start)/v1x;
V_start=a0_Old;  %volt
V_end=0;   %volt
%%
Current_Voltage_Waveform_1stguess=a0_Old+(a1_Old*(Time_1stguess*Old_Multi).^0.5)+a2_Old*Time_1stguess*Old_Multi;
if V_end>V_start
    Total_Time=Time_1stguess(find(Current_Voltage_Waveform_1stguess<V_end,1,'last'));         %Always start at 0
    Time=Time_1stguess(1:find(Current_Voltage_Waveform_1stguess<V_end,1,'last'));
    Current_Voltage_Waveform=Current_Voltage_Waveform_1stguess(1:find(Current_Voltage_Waveform_1stguess<V_end,1,'last'));
else
    Total_Time=Time_1stguess(find(Current_Voltage_Waveform_1stguess>V_end,1,'last'));         %Always start at 0
    Time=Time_1stguess(1:find(Current_Voltage_Waveform_1stguess>V_end,1,'last'));
    Current_Voltage_Waveform=Current_Voltage_Waveform_1stguess(1:find(Current_Voltage_Waveform_1stguess>V_end,1,'last'));
    
end
plot(Time,Current_Voltage_Waveform);

%% 1st order derivative of V (dVdt)
dVdt_Waveform=a1_Old*Time.^(-0.5)+a2_Old;

plot(Time,dVdt_Waveform);


%% Read BinFPC
Folder_Path='D:\171123_PZT calibration related\';
cd(Folder_Path);
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
Bin_Not_Using=[];
FPC(Bin_Not_Using)=NaN;

%% Generating High-resolution FPC waveform

Total_Time_Recorded=length(Bin)*Number_of_Frame_per_Bin/Camera_Frame_Rate;

Time_Bin=Number_of_Frame_per_Bin/Camera_Frame_Rate*(Bin-0.5)*1000;   %ms, move the time tag to the central instead of the end of bin

FPC_hires=interp1(Time_Bin,FPC,Time);

Time_noNaN=Time;
Time_noNaN(isnan(FPC_hires))=[];
FPC_hires_noNaN=FPC_hires;
FPC_hires_noNaN(isnan(FPC_hires))=[];

plot(Time_Bin,FPC,Time_noNaN,FPC_hires_noNaN);

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
Time_New_NoNaN=Time_New_NoNaN*Old_Multi;
plot(Time_New_NoNaN,Voltage_Waveform_New_Equalspacing_noNaN);

%% Fitting (Method 1)
Time_Output_Interp=1:round(Total_Time_Calc*1000/Interp_Target_Multi);
gs = fittype( @(a0,a1, a2, x) a0+a1*x.^0.5+a2*x);
COEFS=fit(Time_New_NoNaN',Voltage_Waveform_New_Equalspacing_noNaN',gs,'StartPoint',[a1_Old a1_Old a2_Old]);%,'StartPoint',[0 0 COEF2.a COEF2.b COEF2.c]);
Voltage_Waveform_Interp=interp1(Time_New_NoNaN/Interp_Target_Multi,Voltage_Waveform_New_Equalspacing_noNaN,Time_Output_Interp,'linear','extrap');
Voltage_Waveform_Interp_No_Ectrap=interp1(Time_New_NoNaN/Interp_Target_Multi,Voltage_Waveform_New_Equalspacing_noNaN,Time_Output_Interp,'linear');
First_NoNAN_Index=find(Voltage_Waveform_Interp_No_Ectrap>0,1,'first');
Voltage_Waveform_Spline=spline(Time_New_NoNaN/Interp_Target_Multi,Voltage_Waveform_New_Equalspacing_noNaN,Time_Output_Interp);

plot(Time_New_NoNaN,Voltage_Waveform_New_Equalspacing_noNaN,Time_New_NoNaN,COEFS.a0+COEFS.a1*Time_New_NoNaN.^0.5+COEFS.a2*Time_New_NoNaN,Time_Output_Interp,Voltage_Waveform_Interp,Time_Output_Interp,Voltage_Waveform_Spline);
legend('Actual','Fit','Interp (Back to 1x)','Spline');

plot(Time_New_NoNaN,Voltage_Waveform_New_Equalspacing_noNaN,Time_New_NoNaN,COEFS.a0+COEFS.a1*Time_New_NoNaN.^0.5+COEFS.a2*Time_New_NoNaN,Time_Output_Interp,Voltage_Waveform_Interp);
legend('Actual','Fit','Interp (Back to 1x)');

plot(Time,COEFS.a0+COEFS.a1*(Time*Old_Multi).^0.5+COEFS.a2*(Time*Old_Multi),Time,a0_Old+a1_Old*(Time*Old_Multi).^0.5+a2_Old*(Time*Old_Multi));
ylim([0 10])


MSE=sum(((COEFS.a0+COEFS.a1*(Time*Old_Multi).^0.5+COEFS.a2*(Time*Old_Multi))-(a0_Old+a1_Old*(Time*Old_Multi).^0.5+a2_Old*(Time*Old_Multi))).^2);
%%
a0_New=10

a1_New=COEFS.a1

a2_New=COEFS.a2

MSE

Last_V=COEFS.a0+COEFS.a1*(max(Time)*Old_Multi).^0.5+COEFS.a2*(max(Time)*Old_Multi)

%%

Voltage_Waveform_Output=Voltage_Waveform_Interp(abs(Voltage_Waveform_Interp-5)<5);


plot(Time_Output_Interp,Voltage_Waveform_Interp);

mkdir([Folder_Path '\Output']);
cd([Folder_Path '\Output']);
dlmwrite(sprintf('%gx_Voltage_Waveform.txt',Interp_Target_Multi),Voltage_Waveform_Output','delimiter','\t','newline','pc','precision', '%.6f');
dlmwrite(sprintf('%gx_Voltage_Waveform_PartialLate23.txt',Interp_Target_Multi),Voltage_Waveform_Output(round(length(Voltage_Waveform_Output)/3:end))','delimiter','\t','newline','pc','precision', '%.6f');
dlmwrite(sprintf('%gx_Voltage_Waveform_PartialLate12.txt',Interp_Target_Multi),Voltage_Waveform_Output(round(length(Voltage_Waveform_Output)/2:end))','delimiter','\t','newline','pc','precision', '%.6f');
dlmwrite(sprintf('%gx_Voltage_Waveform_PartialEarly23.txt',Interp_Target_Multi),Voltage_Waveform_Output(1:round(length(Voltage_Waveform_Output)/3*2))','delimiter','\t','newline','pc','precision', '%.6f');
dlmwrite(sprintf('%gx_Voltage_Waveform_PartialEarly12.txt',Interp_Target_Multi),Voltage_Waveform_Output(1:round(length(Voltage_Waveform_Output)/2))','delimiter','\t','newline','pc','precision', '%.6f');

Voltage_Waveform_Output_revserse=Voltage_Waveform_Output;
Voltage_Waveform_Output_revserse=Voltage_Waveform_Output_revserse(end:-1:1);
dlmwrite(sprintf('%gx_Voltage_Waveform_Reverse.txt',Interp_Target_Multi),Voltage_Waveform_Output_revserse','delimiter','\t','newline','pc','precision', '%.6f');
%% Fit waveform for comparision (3param)
Time_Whole=1:500000;
Voltage_Waveform_Whole=COEFS.a0+COEFS.a1*(Time_Whole*Interp_Target_Multi).^0.5+COEFS.a2*(Time_Whole*Interp_Target_Multi);
Index_Start=find(Voltage_Waveform_Whole<10,1,'first');
Index_Final=find(Voltage_Waveform_Whole>0,1,'last');

Time_Output=Time_Whole(Index_Start:Index_Final);
Voltage_Waveform_Output_Fit=Voltage_Waveform_Whole(Index_Start:Index_Final);

Min_Length=min(length(Voltage_Waveform_Output_Fit),length(Voltage_Waveform_Interp));
MSE_3paramvsInterp=sum((Voltage_Waveform_Output_Fit(1:Min_Length)-Voltage_Waveform_Interp(1:Min_Length)).^2)
VMSE_3paramvsInterp=sum((diff(Voltage_Waveform_Interp(1:Min_Length))-diff(Voltage_Waveform_Output_Fit(1:Min_Length))).^2)

dlmwrite(sprintf('%gx_Voltage_Waveform_3pram_Fit_a0_%g_a1_%g_a2%g.txt',Interp_Target_Multi,COEFS.a0,COEFS.a1,COEFS.a2),Voltage_Waveform_Output_Fit','delimiter','\t','newline','pc','precision', '%.6f');


%% Voltage Velocity Analysis
plot(Time_Output_Interp(2:end),diff(Voltage_Waveform_Interp));

%% Fit Velocity Instead
gs_v = fittype( @(a1, a2, x) 0.5*a1*x.^-(0.5)+a2);
gs_v2 = fittype( @(a1, a2, x) 1.5*a1*x.^(0.5)+a2);

V_Voltage_Waveform_Interp=diff(Voltage_Waveform_Interp);

VMSE_Current=999;
%VMSE_Plot=zeros([1 length(11:(Min_Length-10))]);
for p=50:200
    Range1=[First_NoNAN_Index:p];   %to avoid zero and plane due to extrapolation
    Range2=[(p+1):(Min_Length-1)];%[(p+1):1000];%(Min_Length-1)];

    COEFS_R1=fit(Time_Output_Interp(Range1)',V_Voltage_Waveform_Interp(Range1)',gs_v,'StartPoint',[a1_Old a2_Old]);%,'StartPoint',[0 0 COEF2.a COEF2.b COEF2.c]);
    if If_DifferentModelFor2Seg ==1
        COEFS_R2=fit(Time_Output_Interp(Range2)',V_Voltage_Waveform_Interp(Range2)',gs_v2,'StartPoint',[a1_Old a2_Old]);%,'StartPoint',[0 0 COEF2.a COEF2.b COEF2.c]);
    else
        COEFS_R2=fit(Time_Output_Interp(Range2)',V_Voltage_Waveform_Interp(Range2)',gs_v,'StartPoint',[a1_Old a2_Old]);%,'StartPoint',[0 0 COEF2.a COEF2.b COEF2.c]);
    end
    Curve1=Voltage_Waveform_Interp(1)+COEFS_R1.a1*(Time_Output_Interp).^0.5+COEFS_R1.a2*(Time_Output_Interp);
    LastValueofCurve1=Curve1(p+1);
    if If_DifferentModelFor2Seg ==1
        Curve2=COEFS_R2.a1*(Time_Output_Interp).^(1.5)+COEFS_R2.a2*(Time_Output_Interp);
    else
        Curve2=COEFS_R2.a1*(Time_Output_Interp).^(0.5)+COEFS_R2.a2*(Time_Output_Interp);
    end
    Curve2=Curve2-Curve2(p+1)+LastValueofCurve1;
    Curve_Combined=Curve1;
    Curve_Combined((p+1):end)=Curve2((p+1):end);
    
    VMSE=sum((diff(Voltage_Waveform_Interp(1:Min_Length))-diff(Curve_Combined(1:Min_Length))).^2);
    if VMSE<VMSE_Current
        Curve_Combined_Best=Curve_Combined;
        Curve1_Best=Curve1;
        Curve2_Best=Curve2;
        p_Best=p;
        VMSE_Current=VMSE;
        VMSE_Plot(p)=VMSE;
    end
    disp(p);
end

V_Curve1=diff(Curve1);
V_Curve2=diff(Curve2);
V_Curve_Combined_Best=diff(Curve_Combined_Best);

plot(Time_Output_Interp(2:end),V_Voltage_Waveform_Interp,Time_Output_Interp(2:end),V_Curve1,Time_Output_Interp(2:end),V_Curve2);
legend('Interp','Fit1','Fit2');
ylim([-10E-3 -2E-3]);

plot(Time_Output_Interp(2:end),V_Voltage_Waveform_Interp,Time_Output_Interp(2:end),V_Curve_Combined_Best);
legend('Interp','Combined');
ylim([-10E-3 -2E-3]);
plot(Time_Output_Interp,Voltage_Waveform_Interp,Time_Output_Interp,Curve_Combined);
legend('Interp','2Seg');

V_2Seg_Output=Curve_Combined_Best(abs(Curve_Combined_Best-5)<5);
if If_DifferentModelFor2Seg ==1
    dlmwrite(sprintf('%gx_Voltage_Waveform_2Seg_Seamless_VMSEopt_VelocityFit.txt',Interp_Target_Multi),V_2Seg_Output','delimiter','\t','newline','pc','precision', '%.6f');
else
    dlmwrite(sprintf('%gx_Voltage_Waveform_2Seg_Seamless_VMSEopt_VelocityFit_SameModel.txt',Interp_Target_Multi),V_2Seg_Output','delimiter','\t','newline','pc','precision', '%.6f');
end

%% One Segment Only (New Polynomial)
gs_v_exp=fittype( @(a1, a2, aex, x) (aex)*a1*x.^(aex-1)+a2);
Range2_Best=[(p_Best+1):(Min_Length-1)];%[(p+1):1000];%(Min_Length-1)];

COEFS_EX=fit(Time_Output_Interp(Range2_Best)',V_Voltage_Waveform_Interp(Range2_Best)',gs_v_exp,'StartPoint',[a1_Old a2_Old, 1.5]);%,'StartPoint',[0 0 COEF2.a COEF2.b COEF2.c]);

CurveEX=COEFS_EX.a1*(Time_Output_Interp).^(COEFS_EX.aex)+COEFS_EX.a2*(Time_Output_Interp);
CurveEX=CurveEX-CurveEX(1)+10;
V_CurveEX=diff(CurveEX);
V_2SegEX_Output=CurveEX(abs(CurveEX-5)<5);

plot(Time_Output_Interp(2:end),V_Voltage_Waveform_Interp,Time_Output_Interp(2:end),V_CurveEX);
legend('Interp','Curve EX');
ylim([-10E-3 -2E-3]);
dlmwrite(sprintf('%gx_Voltage_Waveform_CurveEX.txt',Interp_Target_Multi),V_2SegEX_Output','delimiter','\t','newline','pc','precision', '%.6f');
