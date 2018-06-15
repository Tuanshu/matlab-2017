clear all
%% Operations
If_DF=1;
If_FF=0;
If_normalization=1;

If_STFT_1D=1;   % only use the central 1D map for STFT, reducing data szie
Time_Window=100;

Normalization_ROI=[200 550];
Spectrogram_X_ROI=[260 600];
Spectrogram_Y_ROI=[1 8];
RI=1.406;    %Refractive Index
Predicted_Central_Wavelength=0.78;
nBin=1024;
%%

Save_File_Name='170703_Processed Stack.raw';
root_folder_path='C:\TuanShu\170614_PhotonFocus with EMScanner\170620_Numerical Focusing Test\Single Raw 648 8 8128\';
last_folder_name=[1];


If_EM2=1;
If_Auto_Brightness=1;
Auto_Brightness_TH=0.99;
N=4;

If_Numerical_Process=0; %1: for STFT

if If_Numerical_Process == 1
    Temporal_Window=100;
end


If_EffMap=0;%0;
If_Norm=1;

If_Enhace_Deep_Signal=1;


Data_Save_Folder=[root_folder_path '\Processed Data\'];

If_AlterAxialProcess=1;

% Data format related
Row=648;
Colomn=8;
Page=8128;
Byte_Skip=0;
% Processing related
Column_Binning_Factor=1;
Row_Binning_Factor=1;
Lateral_Ave_Factor=1;   %New Param 170413
Axial_Ave_Factor=2*N; %8*8=64, *N for N-point

%Y_Ave_Factor=8;
Product_of_Axial_Decimation_Factor_and_Ave_Factor=round(2);    %should be 64 for 0.28/4, use 8 here
ave_factor=Product_of_Axial_Decimation_Factor_and_Ave_Factor;%[4*4*4];

Phase_Shift=0;

for QQQ=1:length(last_folder_name)
    folder_path=[root_folder_path num2str(last_folder_name(QQQ)) '\'];
    cd(folder_path);
    micron_per_frame=0.2/ave_factor/N;
    Offset_1=0;
    Offset_2=0; 

    Max_Number_of_Frame=3000000;
    if If_EffMap ==1
        Processed_Data_Path=[Data_Save_Folder Save_File_Name];
    elseif If_EffMap ==0
        Processed_Data_Path=[Data_Save_Folder Save_File_Name];
    elseif If_EffMap ==-1
        Processed_Data_Path=[Data_Save_Folder Save_File_Name];
    end
    %%
    file_list=dir(folder_path);    
    for k = length(file_list):-1:1
        if file_list(k).isdir
            file_list(k) = [ ];
            continue
        end
        fname = file_list(k).name;
        if (fname(1) == '.') || strcmp(fname,'Processed Data') ||  strcmp(fname,'Other')
            file_list(k) = [ ];
        end
    end
    %%
    if length(file_list)==1
        fin = fopen(file_list(1).name);
        Image_Temp=fread(fin,[Row,inf],'uint16','b');
        fclose(fin);
        
        %
        Colomn_Total=size(Image_Temp,2);

        Frame=Colomn_Total/Colomn;
        Image_Stack=zeros(Row,Colomn,Frame);
        for r=1:Frame
            Image_Stack(:,:,r)=Image_Temp(:,(1+(r-1)*Colomn):(r*Colomn));
        end
        
        clear Image_Temp
        
    else
        error('Multiple Files in Folder');
    end
    
    Ave_Array=mean(mean(Image_Stack(Normalization_ROI(1):Normalization_ROI(2),:,:),1),2);
   
    if If_normalization ==1
        Image_Stack=Image_Stack./repmat(Ave_Array,[size(Image_Stack,1),size(Image_Stack,2),1])*mean(Ave_Array);
    end
    
    Ave_Image=mean(Image_Stack,3);

    
    if If_DF ==1
        Image_Stack=Image_Stack-repmat(Ave_Image,[1,1,size(Image_Stack,3)]);
    end
    Ave_Array_After_DF=mean(mean(Image_Stack(Normalization_ROI(1):Normalization_ROI(2),:,:),1),2);
    Ave_Image_After_DF=mean(Image_Stack,3);

    if If_FF ==1
        Image_Stack=Image_Stack./repmat(Ave_Image,[1,1,size(Image_Stack,3)])*mean(Ave_Image(:));
    end
    Ave_Array_After_FF=mean(mean(Image_Stack(Normalization_ROI(1):Normalization_ROI(2),:,:),1),2);
    
    if If_STFT_1D
        Y_index=floor(size(Image_Stack,2)/2+0.5);
        Image_Virtual_B_Scan=reshape(Image_Stack(Spectrogram_X_ROI(1):Spectrogram_X_ROI(2),Spectrogram_Y_ROI(1):Spectrogram_Y_ROI(2),:),(Spectrogram_X_ROI(2)-Spectrogram_X_ROI(1)+1)*(Spectrogram_Y_ROI(2)-Spectrogram_Y_ROI(1)+1),[]);
        Image_Spectrogram=zeros(size(Image_Virtual_B_Scan,1),nBin/2+1,floor(size(Image_Virtual_B_Scan,2)/Time_Window));
        for p=1:size(Image_Virtual_B_Scan,1)
            Image_Spectrogram(p,:,:)=spectrogram(Image_Virtual_B_Scan(p,:),Time_Window,0,nBin);  %1st 0 for zero overlapping, nBin for DFT bins, it will outout only half spectrum (512+1)
            disp(p);
        end
    end
end
%%
NNN=656;
imagesc(squeeze(abs(Image_Spectrogram(NNN,:,:))));  %lateral is time


%%

Ave_Spectrogram=squeeze(mean(abs(Image_Spectrogram(Spectrogram_X_ROI(1):Spectrogram_X_ROI(2),:,:)),1));
imagesc(abs(Ave_Spectrogram));  %max frequency = 0.5 * sampling frequency

Axial_Spatial_Sampling_Resolution=Predicted_Central_Wavelength/2/RI/N/ave_factor;                    %micron
Axial_Spatial_Sampling_Resolution_OPD=Predicted_Central_Wavelength/2/N/ave_factor;                    %micron
Axial_Spatial_Sampling_Frequency=1/Axial_Spatial_Sampling_Resolution_OPD;   %Sample/micron-OPD
Half_Axial_Spatial_Sampling_Frequency=0.5*Axial_Spatial_Sampling_Frequency;
Predicted_Central_Frequency=1/(Predicted_Central_Wavelength/2); %based on OPD
Predicted_Central_Bin=512/Half_Axial_Spatial_Sampling_Frequency*Predicted_Central_Frequency;
Max_Spectral_Frequency=3E8*Half_Axial_Spatial_Sampling_Frequency/2*1E6;   %/2 for roundtrip; *1E6 for micron to mm, TiSa ~=3.8E14
Predicted_Spectral_Frequency=3E8/(Predicted_Central_Wavelength)*1E6;
Spectral_Frequency_Array=0:(Max_Spectral_Frequency/512):Max_Spectral_Frequency;
Wavelength_of_Interest_Array=0.5:0.05:1.5; %micron
Spectral_Frequency_of_Interst_Array=3E8./Wavelength_of_Interest_Array*1E6;
Ave_Spectrogram_Based_on_Wavelength=interp2(repmat(1:size(Ave_Spectrogram,2),[size(Ave_Spectrogram,1) 1]),repmat(Spectral_Frequency_Array',[1 size(Ave_Spectrogram,2)]),Ave_Spectrogram,repmat(1:size(Ave_Spectrogram,2),[length(Spectral_Frequency_of_Interst_Array) 1]),repmat(Spectral_Frequency_of_Interst_Array',[1 size(Ave_Spectrogram,2)]));
imagesc(Ave_Spectrogram_Based_on_Wavelength);
set(gca,'YTick',1:2:length(Wavelength_of_Interest_Array),'YTickLabel',downsample(Wavelength_of_Interest_Array,2));
xlabel('Temporal (Axial) Bins');
ylabel('Optical Wavelength (micron, assume RI=1.406)');
%Ave_Spectrogram_Based_on_Wavelength=interp2(repmat(Spectral_Frequency_Array',[1 size(Ave_Spectrogram,2)]),repmat(1:size(Ave_Spectrogram,2),[size(Ave_Spectrogram,1) 1]),Ave_Spectrogram,repmat(Spectral_Frequency_of_Interst_Array',[1 size(Ave_Spectrogram,2)]),repmat(1:size(Ave_Spectrogram,2),[size(Ave_Spectrogram,1) 1]));
%%
fid = fopen(Processed_Data_Path, 'w+');
fwrite(fid, flip(Image_Stack,3), 'double');
fclose(fid);
fclose('all');

