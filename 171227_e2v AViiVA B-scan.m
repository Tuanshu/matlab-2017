clear all
cd('C:\TuanShu\MATLAB\TSLib');
%%
root_folder_path='C:\TuanShu\171229_e2v Test\e2v\';
last_folder_name='1x';
Number=[3];
%% PreAVE calulation
Scanning_Speed=140; %micron/sec
Carrier_Length=0.28;
Frame_Rate=1000000/12.896;
N=4;
Decimation_Factor=1;
PreAVE=round(Frame_Rate/(Scanning_Speed/Carrier_Length)/N/Decimation_Factor);
PostAVE=2;
%%
Auto_Brightness_TH=0.99;
Numerical_Data_Shift=32768; %somehow要shift這樣一個值, 有點像是把signed讀成unshigned
If_normalization=1;
If_DF=1;
If_FF=0;
If_Glass_Search=1;
If_Sharpness=1;
If_Lateral_Sharpness_Only=1;

If_Axial_Bandpass=0;
Wavelength_SP=0.5;  %micron
Wavelength_LP=0.95;  %micron
If_Velocity_Analysis=0;
Ave_PreSpectrogram=2;
Time_Window=200/Ave_PreSpectrogram;
nBin=1024;
Predicted_Central_Wavelength=0.78;
RI=1.4;
Wavelength_of_Interest_Array=0.5:0.01:1.2;%0.6:0.01:1.2; %micron
Spectrogram_X_ROI=[400 600];
%% Some Calculation related to Spectrogram

Axial_Spatial_Sampling_Resolution=Predicted_Central_Wavelength/2/RI/N/PreAVE*Ave_PreSpectrogram;                    %micron
Axial_Spatial_Sampling_Resolution_OPD=Predicted_Central_Wavelength/2/N/PreAVE*Ave_PreSpectrogram;                    %micron
Axial_Spatial_Sampling_Frequency=1/Axial_Spatial_Sampling_Resolution_OPD;   %Sample/micron-OPD
Half_Axial_Spatial_Sampling_Frequency=0.5*Axial_Spatial_Sampling_Frequency; %*0.5 FOR DFT
Predicted_Central_Frequency=1/(Predicted_Central_Wavelength/2); %based on OPD
Predicted_Central_Bin=nBin/2/Half_Axial_Spatial_Sampling_Frequency*Predicted_Central_Frequency;
Max_Spectral_Frequency=3E8*Half_Axial_Spatial_Sampling_Frequency/2*1E6;   %/2 for roundtrip; *1E6 for micron to mm, TiSa ~=3.8E14
Predicted_Spectral_Frequency=3E8/(Predicted_Central_Wavelength)*1E6;

Spectral_Frequency_Array=0:(Max_Spectral_Frequency/(nBin/2)):Max_Spectral_Frequency;
Spectral_Frequency_of_Interst_Array=3E8./Wavelength_of_Interest_Array*1E6;
%%
Multiple_Frame_AVE=16;
If_Lateral_Hilbert=0;
Axial_ROI=[0 350]; %micron
Lateral_ROI=[251 650];  %pixel
%%
Column_Binning_Factor=1;
Row_Binning_Factor=1;
Lateral_Ave_Factor=1;   %New Param 170413
Axial_Ave_Factor=2*N; %8*8=64, *N for N-point
Row=1024;
Colomn=1000;
Phase_Shift=0;
%%
Normalization_ROI=[300 700];
%% Folder Loading
parent_folder_path=[root_folder_path last_folder_name];
folder_list=dir(parent_folder_path);
Original_Lnegth_folder_list=length(folder_list);

for p=1:Original_Lnegth_folder_list
    if folder_list(Original_Lnegth_folder_list+1-p).isdir == 0
        folder_list(Original_Lnegth_folder_list+1-p)=[];
    elseif strcmp(folder_list(Original_Lnegth_folder_list+1-p).name,'.') == 1
        folder_list(Original_Lnegth_folder_list+1-p)=[];
    elseif strcmp(folder_list(Original_Lnegth_folder_list+1-p).name,'..') == 1
        folder_list(Original_Lnegth_folder_list+1-p)=[];
    elseif strcmp(folder_list(Original_Lnegth_folder_list+1-p).name,'Processed Data') == 1
        folder_list(Original_Lnegth_folder_list+1-p)=[];
    elseif strcmp(folder_list(Original_Lnegth_folder_list+1-p).name,'backup') == 1
        folder_list(Original_Lnegth_folder_list+1-p)=[];
    end
end

Data_Save_Folder=[parent_folder_path '\Processed Data\'];

if exist(Data_Save_Folder)==0
    mkdir(Data_Save_Folder);
end

if isempty(Number)==0
    folder_list=folder_list(Number);
end



%%
Sharpness_Array=zeros(length(folder_list),1);
for QQQ=1:length(folder_list)
    folder_path=[parent_folder_path '\' folder_list(QQQ).name '\'];

    file_list_ori=dir(folder_path);
    Original_Length_file_list=length(file_list_ori);

    for p=1:Original_Length_file_list
        [pathstr,name,ext] = fileparts(file_list_ori(Original_Length_file_list+1-p).name);
        if file_list_ori(Original_Length_file_list+1-p).isdir ~= 0
            file_list_ori(Original_Length_file_list+1-p)=[];
        elseif strcmp(file_list_ori(Original_Length_file_list+1-p).name,'.') == 1
            file_list_ori(Original_Length_file_list+1-p)=[];
        elseif strcmp(file_list_ori(Original_Length_file_list+1-p).name,'..') == 1
            file_list_ori(Original_Length_file_list+1-p)=[];
        elseif strcmp(ext,'.PNG')
            file_list_ori(Original_Length_file_list+1-p)=[];
        elseif strcmp(file_list_ori(Original_Length_file_list+1-p).name,'Processed Data') == 1
            file_list_ori(Original_Length_file_list+1-p)=[];
        end
    end
    
    %% To Generate the real File List and Frame_Number_in_File_List Based on Data Size Estimation
    file_list=downsample(file_list_ori,Decimation_Factor,0);
    Frame=length(file_list);
    
       %% File Loading
    Attached_Image=zeros(Row,Colomn*Frame);
    for p=1:Frame
        Attached_Image(:,(1+(p-1)*Colomn):(p*Colomn))=imread([folder_path file_list(p).name],'tiff')'-Numerical_Data_Shift;
        disp(p);
    end
    Ave_Array=mean(Attached_Image(Normalization_ROI(1):Normalization_ROI(2),:),1);
    
     %% Basic Operations 
    if If_normalization ==1
        Attached_Image=Attached_Image./repmat(Ave_Array,[size(Attached_Image,1),1])*mean(Ave_Array);
    end    
    Ave_Image=mean(Attached_Image,2);
    if If_DF ==1
        Attached_Image=Attached_Image-repmat(Ave_Image,[1,size(Attached_Image,2)]);
    end
    
    Ave_Array_After_DF=mean(Attached_Image(Normalization_ROI(1):Normalization_ROI(2),:),1);
    Ave_Image_After_DF=mean(Attached_Image,2);

    if If_FF ==1
        Attached_Image=Attached_Image./repmat(Ave_Image_After_DF,[1,size(Attached_Image,2)])*mean(Ave_Image_After_DF(:));
    end
    Ave_Array_After_FF=mean(Attached_Image(Normalization_ROI(1):Normalization_ROI(2),:),1);
    
    
    %% Axial Bandpass
    if If_Axial_Bandpass == 1
        Axial_SR_Raw=Scanning_Speed/Frame_Rate; %micron
        Optical_Full_BW=3E8/(Axial_SR_Raw*1E-6);
        Frequency_Resolution=Optical_Full_BW/size(Attached_Image,2);
        Frequency_SP=3E8/(Wavelength_LP*(1E-6)/RI/2);
        Frequency_LP=3E8/(Wavelength_SP*(1E-6)/RI/2);
        Index_SP=round(Frequency_SP/Frequency_Resolution);
        Index_LP=round(Frequency_LP/Frequency_Resolution);

        Attached_Image_Axial_FFT=fft(Attached_Image,[],2);
        plot(mean(abs(Attached_Image_Axial_FFT),1))
        Attached_Image_Axial_FFT(:,1:Index_SP)=0;
        Attached_Image_Axial_FFT(:,Index_LP:end)=0;
        Attached_Image=real(ifft(Attached_Image_Axial_FFT,[],2));
        clear Attached_Image_Axial_FFT
    end
    
    %% Velocity Analysis (得要加上DF)
    if If_Velocity_Analysis ==1
        Attached_Image_AVE=TSBinning(Attached_Image,2,Ave_PreSpectrogram);
        Image_Spectrogram=zeros(Spectrogram_X_ROI(2)-Spectrogram_X_ROI(1)+1,nBin/2+1,floor(size(Attached_Image_AVE,2)/Time_Window));
        for p=1:(Spectrogram_X_ROI(2)-Spectrogram_X_ROI(1)+1)
            Image_Spectrogram(p,:,:)=spectrogram(Attached_Image_AVE(Spectrogram_X_ROI(1)+p-1,:),Time_Window,0,nBin);  %1st 0 for zero overlapping, nBin for DFT bins, it will outout only half spectrum (512+1)
            disp(p);
        end
        Ave_Spectrogram=squeeze(mean(abs(Image_Spectrogram),1));
        Ave_Spectrogram_Based_on_Wavelength=interp2(repmat(1:size(Ave_Spectrogram,2),[size(Ave_Spectrogram,1) 1]),repmat(Spectral_Frequency_Array',[1 size(Ave_Spectrogram,2)]),Ave_Spectrogram,repmat(1:size(Ave_Spectrogram,2),[length(Spectral_Frequency_of_Interst_Array) 1]),repmat(Spectral_Frequency_of_Interst_Array',[1 size(Ave_Spectrogram,2)]));
        [Max_Value Max_Index]=max(Ave_Spectrogram_Based_on_Wavelength,[],1);
    end
    
  
%% Lateral Hilbert
    if If_Lateral_Hilbert == 1
        Attached_Image=Attached_Image(Lateral_ROI(2):-1:Lateral_ROI(1),:);
        %%
        
        Rounding_Window=1000;
        Filter=gaussmf(1:3*Rounding_Window,[Rounding_Window round(1.5*Rounding_Window)]);
        Attached_Image_Smooth=conv2(Attached_Image,Filter,'same');
        Attached_Image_Residual=Attached_Image-Attached_Image_Smooth;
        SP=40;
        LP=350;
%%
        Attached_Image_FFT=fft(Attached_Image_Residual,[],1);
        Attached_Image_FFT(1:SP,:)=0;
        Attached_Image_FFT(LP:end,:)=0;
        Attached_Image_Reconst=abs(ifft(Attached_Image_FFT,[],1));
    end
    
    %% PreAVE
   Attached_Image=TSBinning(Attached_Image,2,PreAVE);
   
   %% N-point or Lateral Hilbert
    if If_Lateral_Hilbert == 0
        Temp=0;
        Temp_2=0;
        Length_after_NPoint=floor(size(Attached_Image,2)/N)*N;
        for p=1:N
            Temp=Temp+Attached_Image(:,(N-(p-1)):N:(Length_after_NPoint-(p-1)));
            Temp_2=Temp_2+Attached_Image(:,(N-(p-1)):N:(Length_after_NPoint-(p-1))).^2;
            disp(p);
        end
        Attached_Image=(N*Temp_2-Temp.^2).^0.5*(2^0.5)/N; 
    end

    
    %% PostAVE
    Attached_Image=TSBinning(Attached_Image,2,PostAVE);
    %%
    Ave_Array_After_PostAVE=mean(Attached_Image(Normalization_ROI(1):Normalization_ROI(2),:),1);
    Ave_Array_After_PostAVE_Temp=Ave_Array_After_PostAVE;
    Ave_Array_After_PostAVE_Temp(isnan(Ave_Array_After_PostAVE))=[];
    Ave_Array_After_PostAVE(isnan(Ave_Array_After_PostAVE))=mean(Ave_Array_After_PostAVE_Temp);
    %% Glass Search and  ROI
    if If_Glass_Search == 1
        Glass_TH=10;    %times mean of AVE array
        Glass_Margin=20;    %micron
        if Multiple_Frame_AVE  > 1
            peakfind=1;
            Number_PeakFound=0;
            Last_ROI_LastIndex=1;
            while peakfind
                Glass_Index_Temp=find(Ave_Array_After_PostAVE((Last_ROI_LastIndex+1):end)>Glass_TH*mean(Ave_Array_After_PostAVE),1,'first')+Last_ROI_LastIndex;
                if (Glass_Index_Temp+round((Axial_ROI(2)+Glass_Margin)/(Carrier_Length*PostAVE)))<length(Ave_Array_After_PostAVE)
                    Number_PeakFound=Number_PeakFound+1;
                    Glass_Index(Number_PeakFound)=Glass_Index_Temp;
                    Last_ROI_LastIndex=Glass_Index_Temp+round((Axial_ROI(2)+Glass_Margin)/(Carrier_Length*PostAVE));
                else
                    peakfind=0;
                end
            end
            Axial_ROI_Index_Basic=round((Axial_ROI-Glass_Margin)/(Carrier_Length*PostAVE));
            Number_PeakUsed=min(Multiple_Frame_AVE,Number_PeakFound);
            Temp=zeros((Lateral_ROI(2)-Lateral_ROI(1)+1),(Axial_ROI_Index_Basic(2)-Axial_ROI_Index_Basic(1)+1),Number_PeakUsed);
            for p=1:Number_PeakUsed
                Axial_ROI_Index=Axial_ROI_Index_Basic+Glass_Index(p);
                Temp(:,:,p)=Attached_Image(Lateral_ROI(2):-1:Lateral_ROI(1),Axial_ROI_Index(1):Axial_ROI_Index(2));
                disp(p);
            end
            %%
            Attached_Image=mean(Temp,3);
            %%
        else
            Glass_Index=find(Ave_Array_After_PostAVE>Glass_TH*mean(Ave_Array_After_PostAVE),1,'first'); 
            Axial_ROI_Index=round((Axial_ROI-Glass_Margin)/(Carrier_Length*PostAVE))+Glass_Index;
            Attached_Image=Attached_Image(Lateral_ROI(2):-1:Lateral_ROI(1),Axial_ROI_Index(1):Axial_ROI_Index(2));
        end
    else
        Axial_ROI_Index=round((Axial_ROI)/(Carrier_Length*PostAVE));
        Attached_Image=Attached_Image(Lateral_ROI(2):-1:Lateral_ROI(1),Axial_ROI_Index(1):Axial_ROI_Index(2));
    end
    %Attached_Image=Attached_Image(Lateral_ROI(1):Lateral_ROI(2),Axial_ROI_Index(1):Axial_ROI_Index(2));
    % Lateral Invert
        %%
    imagesc(Attached_Image');
    colormap(gray);
    
            %% Sharpness
    if If_Sharpness ==1
        
        [Gx, Gz]=gradient(Attached_Image);
        if If_Lateral_Sharpness_Only ==1
            Gz=0;
        end
        S=sqrt(Gx.^2+Gz.^2)/mean(mean(Ave_Array_After_FF))*2^12;  %!!!!! 用最原始intensity norm, 有點像NN
        Sharpness_Array(QQQ)=mean(mean(S(:,:)));
        
        
    end
    
    %%  Auto-Brightness
    Auto_Brightness_TH=0.99;
    nbin=250;
    Histogram=hist(Attached_Image(:),nbin);
    Total=sum(Histogram);
    Int_Hist = cumsum(Histogram)/Total;

    C_max=find(Int_Hist>Auto_Brightness_TH,1,'first')/nbin*max(Attached_Image(:));
    C_min=C_max*0.05;

    Attached_Image=(Attached_Image-C_min)/(C_max-C_min);
    Attached_Image(Attached_Image>1)=1;
    Attached_Image(Attached_Image<0)=0;
    %%
    imagesc(Attached_Image');
    colormap(gray);
    imwrite(Attached_Image',[Data_Save_Folder '\' folder_list(QQQ).name '.png'],'png');


end
    %% 