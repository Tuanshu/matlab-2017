clear all
%% Basic Operations
If_normalization=1;
If_DF=1;
If_FF=0;
Normalization_ROI=[400 600];
%% Partial Spectrum
If_Partial_Spectrum=1;
BW=0.04;        %Guassian sig
Offset=0;    %percent
If_Hilbert=1;   %Output the abs (env)
%% Decorrelation along Z
If_Decorrelation=0; %Along Z, ~= along Time
AVE_preDecorr=2;
AVE_postDecorr=16/AVE_preDecorr;
Decorr_Skip=8;
%% Speckle Size Analysis (Axial)
If_Speckle_Size_Full_Frame=1;
If_MassCenter=0;
%% Speckle Size Analysis (Lateral)
If_Speckle_Size_Full_Lateral=0;
%% Phase Map Generation
If_Phase_Map=0;

%% General STFT
Ave_PreSpectrogram=2;
Time_Window=500/Ave_PreSpectrogram;
nBin=2048;
Spectrogram_X_ROI=[1 1024];%[200 800];
Spectrogram_Y_ROI=[1 11];
%% STFT 1D (PZT Speed Analysis)
If_STFT_1D=0;   % only use the central 1D map for STFT, reducing data szie
Wavelength_of_Interest_Array=0.6:0.01:1.2; %micron
RI=1.406;    %Refractive Index
Predicted_Central_Wavelength=0.78;

%% STFT for Dispersion or GVD
If_STFT_GVD=0;

%% N-point Calculation
If_Npoint=0; %Along Z, ~= along Time
N=4;
AVE_preNp=2;
AVE_postNp=1;%16/AVE_preDecorr/N;
Auto_Brightness_TH=0.99;
If_Enhace_Deep_Signal=1;
BI=1;
BF=0.2;
%% If N-point Image Colored by Phase Map
If_Color=0;
Contrast='Kurt';

%% Post N-point Speckle Size Analysis (Axial)
If_PostNpoint_Speckle_Size_Axial=0;
If_MassCenter=0;

%% Lateral Spatial Frequency Analysis (Post N-point or partial spectrum)

If_Lateral_Frequency_Analysis=0;

Lateral_Frequency_ROI=[1 1024];

%% Decorrelation along Y After N-point
If_Decorrelation_Y=0;
AVE_postDecorr_Y=1;
%% Save Raw File
If_Save_Raw=1;
%%
root_folder_path='C:\TuanShu\';
last_folder_name='170718_For Speckle Size Analysis';
Number=[1];

N=4;

% Data format related
Row=1024;
Colomn=11;
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

Number_of_Frame_per_File=2;

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

for QQQ=1:length(folder_list)
    %% File Listing
    folder_path=[parent_folder_path '\' folder_list(QQQ).name '\'];
    file_list_ori=dir(folder_path);
    file_list=reshape(repmat(file_list_ori,[1 Number_of_Frame_per_File])',1,[])';
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
    Frame_Number_in_File_List=reshape(repmat(1:Number_of_Frame_per_File,[length(file_list) 1])',1,[])';
    file_list=downsample(file_list,floor(Product_of_Axial_Decimation_Factor_and_Ave_Factor/ave_factor),Phase_Shift);
    Frame_Number_in_File_List=downsample(Frame_Number_in_File_List,floor(Product_of_Axial_Decimation_Factor_and_Ave_Factor/ave_factor),Phase_Shift);
    Frame=length(file_list);
    
       %% File Loading
    Image_Stack=zeros(Row,Colomn,Frame);
    for p=1:length(file_list)
        fin=fopen([folder_path file_list(p).name]);
        fseek(fin, Byte_Skip+Row*Colomn*2*(Frame_Number_in_File_List(p)-1), 'bof');
        Image_Stack(:,:,p)=fread(fin,[Row,Colomn],'uint16');
        fclose('all');
        disp(p);
    end
        
    Ave_Array=mean(mean(Image_Stack(Normalization_ROI(1):Normalization_ROI(2),:,:),1),2);
       %% Basic Operations 
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
   %% FFT to Frequency Domain
    if (If_Partial_Spectrum==1) || (If_Phase_Map==1)
        LB=0.07;         % ratio
        RB=0.17;
        Image_Stack_FFT=fft(Image_Stack,[],3);
        Image_Stack_FFT(:,:,1:round(size(Image_Stack_FFT,3)*LB))=0;
        Image_Stack_FFT(:,:,round(size(Image_Stack_FFT,3)*RB):end)=0;
    end
    %% Partial Spectrum
    if If_Partial_Spectrum ==1
        Mean_Spectrum=squeeze(mean(mean(abs(Image_Stack_FFT),1),2));    %%only used for center finding
        Mean_Spectrum(Mean_Spectrum<max(Mean_Spectrum)*0.2)=0;          %only >20% peak are considered
        Peak_Index=sum([1:length(Mean_Spectrum)]'.*Mean_Spectrum)/sum(Mean_Spectrum);
        Filter_Function=zeros(1,1,length(Mean_Spectrum));
        Filter_Function(1,1,:)=gaussmf([1:length(Mean_Spectrum)]',[BW*length(Mean_Spectrum) Peak_Index+round(size(Image_Stack_FFT,3)*Offset)]);
        plot(squeeze(mean(mean(abs(Image_Stack_FFT),1),2)).*squeeze(Filter_Function));
        if If_Hilbert == 0
            Image_Stack=2*real(ifft(Image_Stack_FFT.*repmat(Filter_Function,[size(Image_Stack_FFT,1) size(Image_Stack_FFT,2) 1]),[],3));
        else
            Image_Stack=2*abs(ifft(Image_Stack_FFT.*repmat(Filter_Function,[size(Image_Stack_FFT,1) size(Image_Stack_FFT,2) 1]),[],3));
            Image_Stack=flip(Image_Stack,3);
        end
    end
    %% Decorrelation along Z
    if If_Decorrelation == 1
        Temp=0;
        Length_Used_Pre=floor(size(Image_Stack,3)/AVE_preDecorr)*AVE_preDecorr;
        for p=1:AVE_preDecorr
           Temp=Temp+Image_Stack(:,:,(AVE_preDecorr-(p-1)):AVE_preDecorr:(Length_Used_Pre-(p-1)));
        end
        Image_Stack=Temp/AVE_preDecorr;
        
        for p=1:size(Image_Stack,3)
            if p <(size(Image_Stack,3)-Decorr_Skip)
                Image_Stack(:,:,p)=2*(Image_Stack(:,:,p).*Image_Stack(:,:,p+1+Decorr_Skip))./(Image_Stack(:,:,p).^2+Image_Stack(:,:,p+1+Decorr_Skip).^2);
            else
                Image_Stack(:,:,p)=0;
            end
            disp(p);
        end
        Image_Stack=1-Image_Stack;
        
        Temp=0;
        Length_Used_Post=floor(size(Image_Stack,3)/AVE_postDecorr)*AVE_postDecorr;
        for p=1:AVE_postDecorr
           Temp=Temp+Image_Stack(:,:,(AVE_postDecorr-(p-1)):AVE_postDecorr:(Length_Used_Post-(p-1)));
        end
        Image_Stack=Temp/AVE_postDecorr;
        Decorr_Image=squeeze(mean(Image_Stack,2))';
        Decorr_Image=Decorr_Image(size(Decorr_Image,1):-1:1,:);
        Decorr_Image=Decorr_Image/max(Decorr_Image(:));
        Decorr_Image(Decorr_Image>1)=1;
        Decorr_Image(Decorr_Image<0)=0;
        imagesc(Decorr_Image);
        imwrite(Decorr_Image,[Data_Save_Folder sprintf('%s_Decorr_Pre%g_Post%g_Skip%g.png',folder_list(QQQ).name,AVE_preDecorr,AVE_postDecorr,Decorr_Skip)],'png');
    end

       %% Speckle Size Analysis (Axial)
    if If_Speckle_Size_Full_Frame == 1
        %%
        Axial_Window_Size=100;  %Half
        Downsampling_Ratio=100;
        TH=0.2;
        %p=5600;
        Speckle_Contrast_Stack=zeros(size(Image_Stack,1),size(Image_Stack,2),floor(size(Image_Stack,3)/Downsampling_Ratio));
        for p=1:floor(size(Image_Stack,3)/Downsampling_Ratio)
            %%
            p=40;
            Temp_Image=Image_Stack(:,:,max(1,p*Downsampling_Ratio-Axial_Window_Size):min(size(Image_Stack,3),p*Downsampling_Ratio+Axial_Window_Size));
            %Temp_Image=Temp_Image-min(Temp_Image(:));   %for better contrast
            Auto=ifft(conj(fft(Temp_Image,[],3)).*fft(Temp_Image,[],3),[],3); %故意不shift, 因為是symmetric, 可用積分算一半面積
            if If_MassCenter ==0
                Auto=Auto-repmat(min(Auto,[],3),[1 1 size(Auto,3)]);
                Ratio=cumsum(Auto(:,:,1:floor(size(Auto,3)/2)),3)./repmat(sum(Auto(:,:,1:floor(size(Auto,3)/2)),3),[1 1 floor(size(Auto,3)/2)]);       
                [minvalue0 THindex]=min(abs(Ratio-TH),[],3);
                Speckle_Contrast_Stack(:,:,p)=THindex;
            else
                %Auto=Auto-repmat(min(Auto,[],3),[1 1 size(Auto,3)]);
                %Auto=Auto./repmat(max(Auto,[],3),[1 1 size(Auto,3)]);
                Auto=Auto(:,:,1:(Axial_Window_Size+1));
                %Auto=(Auto-repmat(min(Auto,[],3),[1 1 size(Auto,3)]))./(repmat(max(Auto,[],3),[1 1 size(Auto,3)])-repmat(min(Auto,[],3),[1 1 size(Auto,3)]));
                Auto=Auto./repmat(sum(Auto,3),[1 1 size(Auto,3)]);
                Index_Array(1,1,:)=1:(Axial_Window_Size+1);
                %MassCenter=sum(Auto.*repmat(Index_Array,[size(Auto,1) size(Auto,2) 1]),3)./sum(Auto,3);
                MassCenter=mean(Auto.*repmat(Index_Array.^4,[size(Auto,1) size(Auto,2) 1]),3)./mean(Auto.*repmat(Index_Array.^2,[size(Auto,1) size(Auto,2) 1]),3).^2;
                Speckle_Contrast_Stack(:,:,p)=MassCenter;
            end
            disp(p);
            plot(squeeze(Ratio(100,6,:)));
        %%
        end
        
        cmin=5;%6.5;
        cmax=8;%8;
        
        imagesc(squeeze(mean(Speckle_Contrast_Stack(:,:,:),2))');
        caxis([cmin cmax])
        %imagesc(squeeze(mean(Temp_Image(:,:,:),2))');
        %imagesc(fftshift(squeeze(mean(Auto(:,:,:),2))'));
        %plot(fftshift(squeeze(Temp_Image(500,30,:))'));
%%

        %plot(squeeze(mean(mean(Speckle_Contrast_Stack(:,:,:),2),1)));
        
 %% Save Speckle Image
        cmin=5;%6.5;
        cmax=8;%8;
        
        imagesc(squeeze(mean(Speckle_Contrast_Stack(:,:,:),2))');
        caxis([cmin cmax])

        Speckle_Image=squeeze(mean(Speckle_Contrast_Stack(:,:,:),2))';
        Speckle_Image=(Speckle_Image-cmin)./(cmax-cmin);
        Speckle_Image(Speckle_Image>1)=1;
        Speckle_Image(Speckle_Image<0)=0;

        %imwrite(Speckle_Image,[Data_Save_Folder sprintf('%s_Speckle_Image.png',folder_list(QQQ).name)],'png');
        saveas(gcf,[Data_Save_Folder sprintf('%s_Speckle_Image.png',folder_list(QQQ).name)]);

        %% Save Raw of Speckle
        Processed_Data_Path_Speckle=[Data_Save_Folder sprintf('%s_Speckle_DR%g_TH%g.raw',folder_list(QQQ).name,Downsampling_Ratio,TH)];


        fid = fopen(Processed_Data_Path_Speckle, 'w+');
        fwrite(fid, Speckle_Contrast_Stack, 'double');
        fclose(fid);
        fclose('all');
        
%         C_center_mean=mean(C(1,:,1),2)
%         plot(fftshift(mean(C(:,1,1),3),1))
    end
    
    %% Speckle Size Analysis (Lateral)
    if If_Speckle_Size_Full_Lateral == 1
        %%
        Lateral_Window_Size=10;  %Half
        Downsampling_Ratio=10;
        TH=0.5;
        %p=5600;
        Speckle_Contrast_Stack=zeros(floor(size(Image_Stack,1)/Downsampling_Ratio),size(Image_Stack,2),size(Image_Stack,3));
        for p=1:floor(size(Image_Stack,1)/Downsampling_Ratio)
            %%
            %p=16;
            Temp_Image=Image_Stack(max(1,p*Downsampling_Ratio-Lateral_Window_Size):min(size(Image_Stack,1),p*Downsampling_Ratio+Lateral_Window_Size),:,:);
            %Temp_Image=Temp_Image-min(Temp_Image(:));   %for better contrast
            Auto=ifft(conj(fft(Temp_Image,[],1)).*fft(Temp_Image,[],1),[],1); %故意不shift, 因為是symmetric, 可用積分算一半面積
            Auto=Auto-repmat(min(Auto,[],1),[size(Auto,1) 1 1]);
            Ratio=cumsum(Auto(1:floor(size(Auto,1)/2),:,:),1)./repmat(sum(Auto(1:floor(size(Auto,1)/2),:,:),1),[floor(size(Auto,1)/2) 1 1]);       
            [minvalue0 THindex]=min(abs(Ratio-TH),[],1);
            Speckle_Contrast_Stack(p,:,:)=THindex;
            disp(p);
            %plot(squeeze(Temp_Image(500,6,:)));
            plot(squeeze(Auto(:,6,6600)));

        %%
        end
        %imagesc(squeeze(mean(Temp_Image(:,:,:),2))');
        %imagesc(fftshift(squeeze(mean(Auto(:,:,:),2))'));
        %plot(fftshift(squeeze(Temp_Image(500,30,:))'));
        %%
%%

        imagesc(squeeze(mean(Speckle_Contrast_Stack(:,:,:),2))');
        
 %% Save Speckle Image
        cmin=0;%6.5;
        cmax=5;%8;
        
        imagesc(squeeze(mean(Speckle_Contrast_Stack(:,:,:),2))');
        caxis([cmin cmax])

        Speckle_Image=squeeze(mean(Speckle_Contrast_Stack(:,:,:),2))';
        Speckle_Image=(Speckle_Image-cmin)./(cmax-cmin);
        Speckle_Image(Speckle_Image>1)=1;
        Speckle_Image(Speckle_Image<0)=0;

        %imwrite(Speckle_Image,[Data_Save_Folder sprintf('%s_Speckle_Image.png',folder_list(QQQ).name)],'png');
        saveas(gcf,[Data_Save_Folder sprintf('%s_Speckle_Image_Lat.png',folder_list(QQQ).name)]);
        plot(squeeze(mean(mean(Speckle_Contrast_Stack(:,:,:),2),1))');

        %% Save Raw of Speckle
        Processed_Data_Path_Speckle=[Data_Save_Folder sprintf('%s_Speckle_DR%g_TH%g_Lat.raw',folder_list(QQQ).name,Downsampling_Ratio,TH)];


        fid = fopen(Processed_Data_Path_Speckle, 'w+');
        fwrite(fid, Speckle_Contrast_Stack, 'double');
        fclose(fid);
        fclose('all');
        
%         C_center_mean=mean(C(1,:,1),2)
%         plot(fftshift(mean(C(:,1,1),3),1))
    end
    
       %% Phase Map Generation
    if If_Phase_Map ==1
        Image_Stack_Complex=ifft(Image_Stack_FFT,[],3);
        Image_Stack_Phase=unwrap(angle(Image_Stack_Complex),[],3);


        %%
        Ave_PostIPhase=128/8;
        Temp=0;
        Length_Used_PostPhase=floor(size(Image_Stack_Phase,3)/Ave_PostIPhase)*Ave_PostIPhase;
        for p=1:Ave_PostIPhase
           Temp=Temp+Image_Stack_Phase(:,:,(Ave_PostIPhase-(p-1)):Ave_PostIPhase:(Length_Used_PostPhase-(p-1)));
        end
        Image_Stack_Phase_Aved=Temp/Ave_PostIPhase;
        %
        Image_Stack_Phase_Velocity=diff(Image_Stack_Phase_Aved,1,3);
        Image_Stack_Phase_Acc=diff(Image_Stack_Phase_Aved,2,3);
        Image_Stack_Phase_2ndAcc=diff(Image_Stack_Phase_Aved,3,3);

             %%
        NNN=14;
        imagesc(rot90(squeeze((Image_Stack_Phase_Velocity(:,NNN,:)))));
     %%
        STD_Phase_Image=rot90(squeeze(std(Image_Stack_Phase_Velocity,0,2)));
        imagesc(1-STD_Phase_Image);
             
             
             %%

        Ave_LatAve=8/8;
        Temp=0;
        Width_Used_LatAve=floor(size(Image_Stack_Phase_Velocity,1)/Ave_LatAve)*Ave_LatAve;
        for p=1:Ave_LatAve
           Temp=Temp+Image_Stack_Phase_Velocity((Ave_LatAve-(p-1)):Ave_LatAve:(Width_Used_LatAve-(p-1)),:,:);
        end
        Image_Stack_Phase_Velocity_LatAve=Temp/Ave_LatAve;
        Phase_Image=rot90(squeeze(mean(Image_Stack_Phase_Velocity_LatAve,2)));
        %%
        Sigma=10;
        Phase_Image_Latnorm=Phase_Image;
        Phase_Image_Latnorm=Phase_Image./(repmat(mean(Phase_Image,2),[1 size(Phase_Image,2)]));
        Phase_Image_Latnorm_Blur = imgaussfilt(Phase_Image_Latnorm,Sigma);
        Phase_Image_Latnorm_norm=(Phase_Image_Latnorm_Blur-min(Phase_Image_Latnorm_Blur(:)))/(max(Phase_Image_Latnorm_Blur(:))-min(Phase_Image_Latnorm_Blur(:)));
        imagesc(Phase_Image_Latnorm_norm);
        imwrite(Phase_Image_Latnorm_norm,[Data_Save_Folder sprintf('%s_PhaseMap_Lat.png',folder_list(QQQ).name)],'png');

%         %%
%         subplot(1,1,1)
%         NNN=63;
%         imagesc(Image_Stack_Phase(:,:,NNN)-Image_Stack_Phase(:,:,1));
        
        %%
%         XX=48;
%         YY=6;
%         subplot(2,1,1)
%         plot(squeeze(Image_Stack_Phase_Aved(XX,YY,:)));
%         subplot(2,1,2)
%         plot(diff(squeeze(Image_Stack_Phase_Aved(XX,YY,:))));
        %%
    end

       %% Spectrogram for STFT
    if (If_STFT_1D == 1) || If_STFT_GVD
        Image_Virtual_B_Scan=reshape(Image_Stack(Spectrogram_X_ROI(1):Spectrogram_X_ROI(2),Spectrogram_Y_ROI(1):Spectrogram_Y_ROI(2),:),(Spectrogram_X_ROI(2)-Spectrogram_X_ROI(1)+1)*(Spectrogram_Y_ROI(2)-Spectrogram_Y_ROI(1)+1),[]);
        Temp=0;
        Length_Used_PreSpectrogram=floor(size(Image_Virtual_B_Scan,2)/Ave_PreSpectrogram)*Ave_PreSpectrogram;
        for p=1:Ave_PreSpectrogram
           Temp=Temp+Image_Virtual_B_Scan(:,(Ave_PreSpectrogram-(p-1)):Ave_PreSpectrogram:(Length_Used_PreSpectrogram-(p-1)));
        end
        Image_Virtual_B_Scan_Aved=Temp/Ave_PreSpectrogram;
        
        Image_Spectrogram=zeros(size(Image_Virtual_B_Scan_Aved,1),nBin/2+1,floor(size(Image_Virtual_B_Scan_Aved,2)/Time_Window));
        for p=1:size(Image_Virtual_B_Scan_Aved,1)
            Image_Spectrogram(p,:,:)=spectrogram(Image_Virtual_B_Scan_Aved(p,:),Time_Window,0,nBin);  %1st 0 for zero overlapping, nBin for DFT bins, it will outout only half spectrum (512+1)
            disp(p);
        end
    end
       
       %% STFT 1D (PZT Speed Analysis)
    if If_STFT_1D
        Ave_Spectrogram=squeeze(mean(abs(Image_Spectrogram(Spectrogram_X_ROI(1):Spectrogram_X_ROI(2),:,:)),1));
        imagesc(abs(Ave_Spectrogram));  %max frequency = 0.5 * sampling frequency
        caxis([0 500]);
        Axial_Spatial_Sampling_Resolution=Predicted_Central_Wavelength/2/RI/N/ave_factor*Ave_PreSpectrogram;                    %micron
        Axial_Spatial_Sampling_Resolution_OPD=Predicted_Central_Wavelength/2/N/ave_factor*Ave_PreSpectrogram;                    %micron
        Axial_Spatial_Sampling_Frequency=1/Axial_Spatial_Sampling_Resolution_OPD;   %Sample/micron-OPD
        Half_Axial_Spatial_Sampling_Frequency=0.5*Axial_Spatial_Sampling_Frequency; %*0.5 FOR DFT
        Predicted_Central_Frequency=1/(Predicted_Central_Wavelength/2); %based on OPD
        Predicted_Central_Bin=512/Half_Axial_Spatial_Sampling_Frequency*Predicted_Central_Frequency;
        Max_Spectral_Frequency=3E8*Half_Axial_Spatial_Sampling_Frequency/2*1E6;   %/2 for roundtrip; *1E6 for micron to mm, TiSa ~=3.8E14
        Predicted_Spectral_Frequency=3E8/(Predicted_Central_Wavelength)*1E6;
        
        Spectral_Frequency_Array=0:(Max_Spectral_Frequency/512):Max_Spectral_Frequency;
        Spectral_Frequency_of_Interst_Array=3E8./Wavelength_of_Interest_Array*1E6;
        Ave_Spectrogram_Based_on_Wavelength=interp2(repmat(1:size(Ave_Spectrogram,2),[size(Ave_Spectrogram,1) 1]),repmat(Spectral_Frequency_Array',[1 size(Ave_Spectrogram,2)]),Ave_Spectrogram,repmat(1:size(Ave_Spectrogram,2),[length(Spectral_Frequency_of_Interst_Array) 1]),repmat(Spectral_Frequency_of_Interst_Array',[1 size(Ave_Spectrogram,2)]));
        [Max_Value Max_Index]=max(Ave_Spectrogram_Based_on_Wavelength,[],1);
        Ave_Spectrogram_Based_on_Wavelength_Norm=Ave_Spectrogram_Based_on_Wavelength./repmat(Max_Value,[size(Ave_Spectrogram_Based_on_Wavelength,1) 1]);
        
        Bin_FPC=[(1:size(Ave_Spectrogram_Based_on_Wavelength_Norm,2))' (Wavelength_of_Interest_Array(Max_Index)/Predicted_Central_Wavelength*N*ave_factor)'];
        
        dlmwrite([Data_Save_Folder sprintf('%s_Bin_FPC.txt',folder_list(QQQ).name)],Bin_FPC,'delimiter','\t','newline','pc','precision', '%.6f');

        imagesc(Ave_Spectrogram_Based_on_Wavelength_Norm);
        colormap(parula);
        hold on
        plot(1:length(Max_Index),Max_Index,'ro');
        hold off
        set(gca,'YTick',1:8:length(Wavelength_of_Interest_Array),'YTickLabel',downsample(Wavelength_of_Interest_Array,8),'FontSize', 14);
        xlabel('Temporal (Axial) Bins');
        ylabel('Optical Wavelength (micron, assume RI=1.406)');
        pbaspect([2 1 1]);

        saveas(gcf,[Data_Save_Folder sprintf('%s_TFA.png',folder_list(QQQ).name)]);
        set(gca,'YTick',(1:8:length(Wavelength_of_Interest_Array)),'YTickLabel',round(downsample(Wavelength_of_Interest_Array/Predicted_Central_Wavelength*N*ave_factor,8),1),'FontSize', 14);
        xlabel('Temporal (Axial) Bins');
        ylabel('Frame per Carrier');
        pbaspect([2 1 1]);
        saveas(gcf,[Data_Save_Folder sprintf('%s_FPC.png',folder_list(QQQ).name)]);
    end
       %% STFT for Dispersion GVD
    if If_STFT_GVD==1
        %%
        Image_Spectrogram_norm=(abs(Image_Spectrogram)-repmat(min(abs(Image_Spectrogram),[],2),[1 size(Image_Spectrogram,2) 1]))./(repmat(max(abs(Image_Spectrogram),[],2),[1 size(Image_Spectrogram,2) 1])-repmat(min(abs(Image_Spectrogram),[],2),[1 size(Image_Spectrogram,2) 1]));
         Image_Spectrogram_MassCenter=sum(abs(Image_Spectrogram_norm).*repmat([1:size(Image_Spectrogram_norm,2)],[size(Image_Spectrogram_norm,1) 1 size(Image_Spectrogram_norm,3)]),2)./sum(abs(Image_Spectrogram_norm),2);

        
        %% Find a way to find the dispersed BW
        NNN=1981;
        imagesc(abs(squeeze(Image_Spectrogram_norm(NNN,:,:))));
        hold on
        plot(1:size(Image_Spectrogram_norm,3),squeeze(Image_Spectrogram_MassCenter(NNN,1,:)),'ro');
        hold off
        
        % try kurtosis
        Index_Array(1,:,1)=1:size(Image_Spectrogram_norm,2);
        Image_Index=repmat(Index_Array,[size(Image_Spectrogram_norm,1) 1 size(Image_Spectrogram_norm,3)]);
        Image_Spectrogram_Kurt=sum((Image_Spectrogram_norm.*(Image_Index-repmat(Image_Spectrogram_MassCenter,[1 size(Image_Spectrogram_norm,2) 1]))).^4,2)./sum((Image_Spectrogram_norm.*(Image_Index-repmat(Image_Spectrogram_MassCenter,[1 size(Image_Spectrogram_norm,2) 1])).^2).^2,2);
        imagesc(squeeze(Image_Spectrogram_Kurt));
        Image_Kurt_Bscan=flip(squeeze(mean(reshape(Image_Spectrogram_Kurt,(Spectrogram_X_ROI(2)-Spectrogram_X_ROI(1)+1),[],size(Image_Spectrogram_Kurt,3)),2))',1);
        imagesc(squeeze(Image_Kurt_Bscan));

        %%
        %imagesc(squeeze(abs(Image_Spectrogram(50,:,:))));  %max frequency = 0.5 * sampling frequency
        Image_Spectrogram_Phase=unwrap(angle(Image_Spectrogram),[],2);
        Image_Spectrogram_1stDerivative=diff(Image_Spectrogram_Phase,1,2);
        Image_Spectrogram_2ndDerivative=diff(Image_Spectrogram_Phase,2,2);
        %Image_Spectrogram_3rdDerivative=diff(Image_Spectrogram_Phase,3,2);
        %imagesc(squeeze(real(Image_Spectrogram_2ndDerivative(50,:,:))));  %max frequency = 0.5 * sampling frequency
        Image_1stDerivative=squeeze(sum(abs(Image_Spectrogram(:,1:(size(Image_Spectrogram,2)-1),:)).*Image_Spectrogram_1stDerivative,2)./sum(abs(Image_Spectrogram(:,1:(size(Image_Spectrogram,2)-1),:)),2));
        Image_2ndDerivative=squeeze(sum(abs(Image_Spectrogram(:,1:(size(Image_Spectrogram,2)-2),:)).*Image_Spectrogram_2ndDerivative,2)./sum(abs(Image_Spectrogram(:,1:(size(Image_Spectrogram,2)-2),:)),2));
        %Image_3rdDerivative=squeeze(sum(abs(Image_Spectrogram(:,1:(size(Image_Spectrogram,2)-3),:)).*Image_Spectrogram_3rdDerivative,2)./sum(abs(Image_Spectrogram(:,1:(size(Image_Spectrogram,2)-3),:)),2));

        Image_1stDerivative_Bscan=flip(squeeze(mean(reshape(Image_1stDerivative,(Spectrogram_X_ROI(2)-Spectrogram_X_ROI(1)+1),[],size(Image_1stDerivative,2)),2))',1);
        Image_2ndDerivative_Bscan=flip(squeeze(mean(reshape(Image_2ndDerivative,(Spectrogram_X_ROI(2)-Spectrogram_X_ROI(1)+1),[],size(Image_2ndDerivative,2)),2))',1);
        %Image_3rdDerivative_Bscan=flip(squeeze(mean(reshape(Image_3rdDerivative,(Spectrogram_X_ROI(2)-Spectrogram_X_ROI(1)+1),[],size(Image_2ndDerivative,2)),2))',1);
        Image_1stDerivative_Volume=reshape(Image_1stDerivative,(Spectrogram_X_ROI(2)-Spectrogram_X_ROI(1)+1),[],size(Image_1stDerivative,2));
        Image_2ndDerivative_Volume=reshape(Image_2ndDerivative,(Spectrogram_X_ROI(2)-Spectrogram_X_ROI(1)+1),[],size(Image_2ndDerivative,2));
%%
%         %%
%         Image_Spectrogram_MassCenter_Bscan=flip(squeeze(mean(reshape(Image_Spectrogram_MassCenter,(Spectrogram_X_ROI(2)-Spectrogram_X_ROI(1)+1),[],size(Image_Spectrogram_MassCenter,3)),2))',1);
%         
%         Image_Spectrogram_MassCenter_Bscan_LatNorm=Image_Spectrogram_MassCenter_Bscan./(repmat(mean(Image_Spectrogram_MassCenter_Bscan,2),[1 size(Image_Spectrogram_MassCenter_Bscan,2)]));
%         Image_Spectrogram_MassCenter_Bscan_LatNorm_norm=(Image_Spectrogram_MassCenter_Bscan_LatNorm-min(Image_Spectrogram_MassCenter_Bscan_LatNorm(:)))/(max(Image_Spectrogram_MassCenter_Bscan_LatNorm(:))-min(Image_Spectrogram_MassCenter_Bscan_LatNorm(:)));
%         imagesc(Image_Spectrogram_MassCenter_Bscan_LatNorm_norm);
        
        
        subplot(2,1,1)
        imagesc((Image_1stDerivative_Bscan))
        %%imagesc(Image_Spectrogram_MassCenter_Bscan)


        colormap(gray)
        %%
        %caxis([0.5 1.5]);
        subplot(2,1,2)
        imagesc(abs(Image_2ndDerivative_Bscan))
        colormap(gray)
        %subplot(2,1,3)
        %imagesc((Image_3rdDerivative_Bscan))
    end
    %% N-point
    if If_Npoint == 1
        Temp=0;
        Length_Used_Pre=floor(size(Image_Stack,3)/AVE_preNp)*AVE_preNp;
        for p=1:AVE_preNp
           Temp=Temp+Image_Stack(:,:,(AVE_preNp-(p-1)):AVE_preNp:(Length_Used_Pre-(p-1)));
        end
        Image_Stack=Temp/AVE_preNp;
        
        Temp=0;
        Temp_2=0;
        Length_Used_Pre=floor(size(Image_Stack,3)/N)*N;
        for p=1:N
            Temp=Temp+Image_Stack(:,:,(N-(p-1)):N:(Length_Used_Pre-(p-1)));
            Temp_2=Temp_2+Image_Stack(:,:,(N-(p-1)):N:(Length_Used_Pre-(p-1))).^2;
            disp(p);
        end
        Image_Stack=(N*Temp_2-Temp.^2).^0.5*(2^0.5)/N; 
        
        Temp=0;
        Length_Used_Post=floor(size(Image_Stack,3)/AVE_postNp)*AVE_postNp;
        for p=1:AVE_postNp
           Temp=Temp+Image_Stack(:,:,(AVE_postNp-(p-1)):AVE_postNp:(Length_Used_Post-(p-1)));
        end
        Image_Stack=Temp/AVE_postNp;   
        Image_Stack=Image_Stack(:,:,size(Image_Stack,3):-1:1);

        
%%
        if If_Enhace_Deep_Signal == 1

            Linear_Axial_Profile=[BI:-1*(BI-BF)/(size(Image_Stack,3)-1):BF]';
    %         imagesc(Reduced_Image)
    %         caxis([C_min C_max])

    %         Profile_Smoothing_Window=10;
    % 
    %         Image_Axial_Profile=mean(Reduced_Image,2);
    % 
    %         Image_Axial_Profile_Smooth=conv(Image_Axial_Profile,ones(Profile_Smoothing_Window,1)/Profile_Smoothing_Window,'same');
    %         plot(Image_Axial_Profile_Smooth)
            Enhace_Coef(1,1,:)=mean(Linear_Axial_Profile)./Linear_Axial_Profile;

            Image_Stack=Image_Stack.*repmat(Enhace_Coef,[size(Image_Stack,1) size(Image_Stack,2) 1]);

        end   
        
        %%
        Npoint_Image=squeeze(mean(Image_Stack,2))';

     
 
        nbin=250;
        Histogram=hist(Npoint_Image(:),nbin);
        Total=sum(Histogram);
        Int_Hist = cumsum(Histogram)/Total;

        C_max=find(Int_Hist>Auto_Brightness_TH,1,'first')/nbin*max(Npoint_Image(:));
        C_min=C_max*0.05;

        Npoint_Image=(Npoint_Image-C_min)/(C_max-C_min);
        Npoint_Image(Npoint_Image>1)=1;
        Npoint_Image(Npoint_Image<0)=0;
        imagesc(Npoint_Image);
        colormap(gray);
        
%         %% Lat norm test
%         Sigma=10;
%         Npoint_Image_LatNorm=Npoint_Image./(repmat(mean(Npoint_Image,2),[1 size(Npoint_Image,2)]));
%         Npoint_Image_LatNorm_Blur = imgaussfilt(Npoint_Image_LatNorm,Sigma);
%         Npoint_Image_LatNorm_norm=(Npoint_Image_LatNorm_Blur-min(Npoint_Image_LatNorm_Blur(:)))/(max(Npoint_Image_LatNorm_Blur(:))-min(Npoint_Image_LatNorm_Blur(:)));
%         imagesc(Npoint_Image_LatNorm_norm);
%         imwrite(Npoint_Image_LatNorm_norm,[Data_Save_Folder sprintf('%s_Npoint_Intensity Latnorm.png',folder_list(QQQ).name)],'png');
    %%
%     Phase_Image_Latnorm_norm_divbyNpoint=(Phase_Image_Latnorm_norm)-(Npoint_Image_LatNorm_norm(1:(size(Npoint_Image_LatNorm_norm)-1),:));
%                 imagesc(Phase_Image_Latnorm_norm);
%     caxis([0 1]);
        %%
        if If_Color == 1

            if strcmp(Contrast,'PhaseMap')
                Color_Map_Ori=Phase_Image_Latnorm_norm;
            elseif strcmp(Contrast,'1stDerivative')
                Color_Map_Ori=Image_1stDerivative_Bscan;      
                
            elseif strcmp(Contrast,'PhaseMap_New')
                Color_Map_Ori=Phase_Image_Latnorm_norm_divbyNpoint;  Kurt
                
            elseif strcmp(Contrast,'Kurt')
                Color_Map_Ori=Image_Kurt_Bscan; 
                
            end
            Npoint_Height_Used=floor(size(Npoint_Image,1)/size(Color_Map_Ori,1))*size(Color_Map_Ori,1);
            Npoint_Width_Used=floor(size(Npoint_Image,2)/size(Color_Map_Ori,2))*size(Color_Map_Ori,2);
            Intensity_Map=Npoint_Image(1:Npoint_Height_Used,1:Npoint_Width_Used);
            Color_Map_2D=interp2(repmat([0:(size(Color_Map_Ori,2)+1)/size(Color_Map_Ori,2):size(Color_Map_Ori,2)],[size(Color_Map_Ori,1) 1]),repmat([0:(size(Color_Map_Ori,1)+1)/size(Color_Map_Ori,1):size(Color_Map_Ori,1)]',[1 size(Color_Map_Ori,2)]),Color_Map_Ori,repmat([1:size(Intensity_Map,2)],[size(Intensity_Map,1) 1])/size(Intensity_Map,2)*size(Color_Map_Ori,2),repmat([1:size(Intensity_Map,1)]',[1 size(Intensity_Map,2)])/size(Intensity_Map,1)*size(Color_Map_Ori,1),'linear',mean(Color_Map_Ori(:)));

            imagesc(Color_Map_2D);
            
            Color_Map_2D_norm=(Color_Map_2D-min(Color_Map_2D(:)))/(max(Color_Map_2D(:))-min(Color_Map_2D(:)));
            %%
            Color_1=1*[1.2 0 0];
            Color_2=1*[0 0.4 0.2];
            
            Color_Map_3D=zeros(size(Color_Map_2D_norm,1),size(Color_Map_2D_norm,2),3);
            for p=1:3
                Color_Map_3D(:,:,p)=Color_Map_2D_norm*Color_1(p)+(1-Color_Map_2D_norm)*Color_2(p);
            end
            Color_Map_3D=Color_Map_3D./repmat(mean(Color_Map_3D,3),[1 1 3]);
            %imagesc(Intensity_Map);
            imagesc(Color_Map_2D_norm);
            imagesc(Color_Map_3D);
            %
            Npoint_Colored=Color_Map_3D.*repmat(Intensity_Map,[1 1 3]);
            
            imagesc(Npoint_Colored);
%%
            imwrite(Npoint_Colored,[Data_Save_Folder sprintf('%s_Npoint_Color.png',folder_list(QQQ).name)],'png');
            imwrite(uint16(Npoint_Colored*2^16),[Data_Save_Folder sprintf('%s_Npoint_Color.tif',folder_list(QQQ).name)],'tiff');
          
        elseif If_Partial_Spectrum == 0
            imwrite(Npoint_Image,[Data_Save_Folder sprintf('%s_Npoint.png',folder_list(QQQ).name)],'png');
            imwrite(uint16(Npoint_Image*2^16),[Data_Save_Folder sprintf('%s_Npoint.tif',folder_list(QQQ).name)],'tiff');

        else
            imwrite(Npoint_Image,[Data_Save_Folder sprintf('%s_Npoint_Filtered_BW%g_Offset%g.png',folder_list(QQQ).name,BW,Offset)],'png');
            imwrite(uint16(Npoint_Image*2^16),[Data_Save_Folder sprintf('%s_Npoint_Filtered_BW%g_Offset%g.tif',folder_list(QQQ).name,BW,Offset)],'tiff');
        end

    end
    
    %% Post N-point Speckle Axial
    if If_PostNpoint_Speckle_Size_Axial == 1
        %%
        clear Index_Array
        Axial_Window_Size=25;  %Half
        Downsampling_Ratio=1;
        TH=0.3;
        %p=5600;
        Speckle_Contrast_Stack=zeros(size(Image_Stack,1),size(Image_Stack,2),floor(size(Image_Stack,3)/Downsampling_Ratio));
        for p=1:floor(size(Image_Stack,3)/Downsampling_Ratio)
            %%
            %p=620;
            Temp_Image=Image_Stack(:,:,max(1,p*Downsampling_Ratio-Axial_Window_Size):min(size(Image_Stack,3),p*Downsampling_Ratio+Axial_Window_Size));
            %Temp_Image=Temp_Image-min(Temp_Image(:));   %for better contrast
            Auto=ifft(conj(fft(Temp_Image,[],3)).*fft(Temp_Image,[],3),[],3); %故意不shift, 因為是symmetric, 可用積分算一半面積
            if If_MassCenter ==0
                Auto=Auto-repmat(min(Auto,[],3),[1 1 size(Auto,3)]);
                Ratio=cumsum(Auto(:,:,1:floor(size(Auto,3)/2)),3)./repmat(sum(Auto(:,:,1:floor(size(Auto,3)/2)),3),[1 1 floor(size(Auto,3)/2)]);       
                [minvalue0 THindex]=min(abs(Ratio-TH),[],3);
                Speckle_Contrast_Stack(:,:,p)=THindex;
            else
                %Auto=Auto-repmat(min(Auto,[],3),[1 1 size(Auto,3)]);
                %Auto=Auto./repmat(max(Auto,[],3),[1 1 size(Auto,3)]);
                Auto=Auto(:,:,1:(Axial_Window_Size+1));
                %Auto=(Auto-repmat(min(Auto,[],3),[1 1 size(Auto,3)]))./(repmat(max(Auto,[],3),[1 1 size(Auto,3)])-repmat(min(Auto,[],3),[1 1 size(Auto,3)]));
                %Auto=Auto./repmat(sum(Auto,3),[1 1 size(Auto,3)]);
                Auto=(Auto-repmat(min(Auto,[],3),[1 1 size(Auto,3)]))./(repmat(max(Auto,[],3),[1 1 size(Auto,3)])-repmat(min(Auto,[],3),[1 1 size(Auto,3)]));
                Auto(Auto<TH)=0;
                Index_Array(1,1,:)=1:(Axial_Window_Size+1);
                %MassCenter=sum(Auto.*repmat(Index_Array,[size(Auto,1) size(Auto,2) 1]),3)./sum(Auto,3);
                MassCenter=mean(Auto.*repmat(Index_Array.^4,[size(Auto,1) size(Auto,2) 1]),3)./mean(Auto.*repmat(Index_Array.^2,[size(Auto,1) size(Auto,2) 1]),3).^2;
                Speckle_Contrast_Stack(:,:,p)=MassCenter;
            end
            disp(p);
            %plot(squeeze(Auto(500,6,:)));
        %%
        end
        %imagesc(squeeze(mean(Temp_Image(:,:,:),2))');
        %imagesc(fftshift(squeeze(mean(Auto(:,:,:),2))'));
        %plot(fftshift(squeeze(Temp_Image(500,30,:))'));
%%

        %plot(squeeze(mean(mean(Speckle_Contrast_Stack(:,:,:),2),1)));
        
 %% Save Speckle Image
        cmin=46;%6.5;
        cmax=50;%8;
        
        imagesc(-squeeze(mean(Speckle_Contrast_Stack(:,:,:),2))');
        caxis([cmin cmax])

        Speckle_Image=squeeze(mean(Speckle_Contrast_Stack(:,:,:),2))';
        Speckle_Image=(Speckle_Image-cmin)./(cmax-cmin);
        Speckle_Image(Speckle_Image>1)=1;
        Speckle_Image(Speckle_Image<0)=0;

        %imwrite(Speckle_Image,[Data_Save_Folder sprintf('%s_Speckle_Image.png',folder_list(QQQ).name)],'png');
        saveas(gcf,[Data_Save_Folder sprintf('%s_Speckle_Image.png',folder_list(QQQ).name)]);

        %% Save Raw of Speckle
        Processed_Data_Path_Speckle=[Data_Save_Folder sprintf('%s_Speckle_DR%g_TH%g.raw',folder_list(QQQ).name,Downsampling_Ratio,TH)];


        fid = fopen(Processed_Data_Path_Speckle, 'w+');
        fwrite(fid, Speckle_Contrast_Stack, 'double');
        fclose(fid);
        fclose('all');
        
%         C_center_mean=mean(C(1,:,1),2)
%         plot(fftshift(mean(C(:,1,1),3),1))
    end
    
%% Lateral Spatial Frequency Analysis (Post N-point or partial spectrum)
if If_Lateral_Frequency_Analysis==1
    Lateral_Frequency_ROI=[200 800];
    Image_Stack_DepthNorm=(Image_Stack(Lateral_Frequency_ROI(1):Lateral_Frequency_ROI(2),:,:)-repmat(min(Image_Stack(Lateral_Frequency_ROI(1):Lateral_Frequency_ROI(2),:,:),[],1),[Lateral_Frequency_ROI(2)-Lateral_Frequency_ROI(1)+1 1 1]))./(repmat(max(Image_Stack(Lateral_Frequency_ROI(1):Lateral_Frequency_ROI(2),:,:),[],1),[Lateral_Frequency_ROI(2)-Lateral_Frequency_ROI(1)+1 1 1])-repmat(min(Image_Stack(Lateral_Frequency_ROI(1):Lateral_Frequency_ROI(2),:,:),[],1),[Lateral_Frequency_ROI(2)-Lateral_Frequency_ROI(1)+1 1 1]));
%%
    imagesc(squeeze(Image_Stack_DepthNorm(:,6,:))');
%%
    Image_Stack_FFTX=fft(Image_Stack_DepthNorm,[],1);
    
    %%
    FFTX_vs_Depth=squeeze(mean(abs(Image_Stack_FFTX(:,:,:)),2));
    FFTX_vs_Depth=FFTX_vs_Depth-repmat(min(FFTX_vs_Depth,[],1),[size(FFTX_vs_Depth,1) 1]);
    %%
    
    plot(FFTX_vs_Depth(:,6504))
    ylim([0 20])
    xlim([0 Lateral_Frequency_ROI(2)-Lateral_Frequency_ROI(1)])

    %%
    Low_Freq_Band=[50 150];
    High_Freq_Band=[150 250];
    Low_Freq_vs_Depth=sum(FFTX_vs_Depth(Low_Freq_Band(1):Low_Freq_Band(2),:),1);
    High_Freq_vs_Depth=sum(FFTX_vs_Depth(High_Freq_Band(1):High_Freq_Band(2),:),1);
    LH_Ratio_vs_Depth=High_Freq_vs_Depth./Low_Freq_vs_Depth;
    plot(LH_Ratio_vs_Depth)
end

        %% Decorrelation along Y after N-oint
    if If_Decorrelation_Y == 1
        
        for p=1:size(Image_Stack,2)
            if p <(size(Image_Stack,2))
                Image_Stack(:,p,:)=2*(Image_Stack(:,p,:).*Image_Stack(:,p+1,:))./(Image_Stack(:,p,:).^2+Image_Stack(:,p+1,:).^2);
            else
                Image_Stack(:,p,:)=0;
            end
            disp(p);
        end
        Image_Stack=1-Image_Stack;
        
        Temp=0;
        Length_Used_Post=floor(size(Image_Stack,3)/AVE_postDecorr_Y)*AVE_postDecorr_Y;
        for p=1:AVE_postDecorr_Y
           Temp=Temp+Image_Stack(:,:,(AVE_postDecorr_Y-(p-1)):AVE_postDecorr_Y:(Length_Used_Post-(p-1)));
        end
        Image_Stack=Temp/AVE_postDecorr_Y;
        
        Decorr_Image=squeeze(mean(Image_Stack,2))';
        Decorr_Image=Decorr_Image(size(Decorr_Image,1):-1:1,:);
        Decorr_Image=Decorr_Image/max(Decorr_Image(:));
        Decorr_Image(Decorr_Image>1)=1;
        Decorr_Image(Decorr_Image<0)=0;
        imagesc(Decorr_Image);
        imwrite(Decorr_Image,[Data_Save_Folder sprintf('%s_Decorr_Y.png',folder_list(QQQ).name)],'png');
    end
    
    %% Save Raw File
    if If_Save_Raw == 1
        if If_Decorrelation == 0
            Processed_Data_Path=[Data_Save_Folder sprintf('%s.raw',folder_list(QQQ).name)];
        else
            Processed_Data_Path=[Data_Save_Folder sprintf('%s_Decorr_Pre%g_Post%g.raw',folder_list(QQQ).name,AVE_preDecorr,AVE_postDecorr)];            
        end

        fid = fopen(Processed_Data_Path, 'w+');
        fwrite(fid, Image_Stack, 'double');
        fclose(fid);
        fclose('all');
    end
end



