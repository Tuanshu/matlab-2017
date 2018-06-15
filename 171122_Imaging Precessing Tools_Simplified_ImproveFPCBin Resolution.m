clear all
%%
Task='PZT';    %CCDZ or PZT
If_New_Folder_Structure=1;
%% Basic Operations
If_normalization=1;
If_DF=1;
If_FF=0;
Normalization_ROI=[400 600];
%% Partial Spectrum
if strcmp(Task,'CCDZ')
    If_Partial_Spectrum=1;
    BW=0.32;        %Guassian sig
    Offset=0;    %percent
    If_Hilbert=1;   %Output the abs (env)
elseif strcmp(Task,'PZT')
    If_Partial_Spectrum=0;
    BW=0.32;        %Guassian sig
    Offset=0;    %percent
    If_Hilbert=0;   %Output the abs (env)
end
%% Image Sharpness Estimation
if strcmp(Task,'CCDZ')
    If_Sharpness_Est=1;
    If_Lateral_Sharpness_Only=1;
    AVE_postHilbert=4*2*2;
    Axial_Resolution=0.56;  %
    Axial_ROI=[100 225];%[100 350];    %micron
    elseif strcmp(Task,'PZT')
    If_Sharpness_Est=0;
    If_Lateral_Sharpness_Only=1;
    AVE_postHilbert=4*2*2;
    Axial_Resolution=0.56;  %
    Axial_ROI=[100 350];    %micron
end


%% General STFT
Ave_PreSpectrogram=2;
Time_Window=200/Ave_PreSpectrogram;
nBin=1024;
Spectrogram_X_ROI=[1 1024];%[200 800];
Spectrogram_Y_ROI=[1 11];
%% STFT 1D (PZT Speed Analysis)
if strcmp(Task,'CCDZ')
    If_STFT_1D=0;   % only use the central 1D map for STFT, reducing data szie
    Wavelength_of_Interest_Array=0.6:0.01:1.2; %micron
    RI=1.406;    %Refractive Index
    Predicted_Central_Wavelength=0.78;
elseif strcmp(Task,'PZT')
    If_STFT_1D=1;   % only use the central 1D map for STFT, reducing data szie
    Wavelength_of_Interest_Array=0.5:0.01:1.2;%0.6:0.01:1.2; %micron
    RI=1.406;    %Refractive Index
    Predicted_Central_Wavelength=0.78;
end



%% N-point Calculation
If_Npoint=0; %Along Z, ~= along Time
N=4;
AVE_preNp=2;
AVE_postNp=2;%16/AVE_preDecorr/N;
Auto_Brightness_TH=0.99;
If_Enhace_Deep_Signal=1;

BI=-0.1;
BF=1;

%% Save Raw File
If_Save_Raw=0;
%%
root_folder_path='F:\AMO\171122_PZT calibration Test\';
last_folder_name='CurveEX Final';
Number=[15:16];

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
ave_factor=2;%[4*4*4];

Phase_Shift=0;

Number_of_Frame_per_File=2;


% For Last File Number Estimation
Bit_per_Pixel=16;
Frame_Size=Row*Colomn*Bit_per_Pixel/8;  %Byte


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
    if If_New_Folder_Structure ==1
        folder_path=[parent_folder_path '\' folder_list(QQQ).name '\Raw Image\'];
    else
        folder_path=[parent_folder_path '\' folder_list(QQQ).name '\'];
    end

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
    file_list=[];
    Frame_Number_in_File_List=[];
    for p=1:length(file_list_ori)
        Frame_Number_Calculated_Based_on_Size=file_list_ori(p).bytes/Frame_Size;
        file_list=[file_list;repmat(file_list_ori(p),[Frame_Number_Calculated_Based_on_Size 1])]; 
        Frame_Number_in_File_List=[Frame_Number_in_File_List;[1:Frame_Number_Calculated_Based_on_Size]']; 
    end
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
    elseif If_normalization ==2
        Image_Stack=Image_Stack./Image_Stack;
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
    if (If_Partial_Spectrum==1)
         LB=0.02;         % ratio
         RB=0.5;
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
        Image_Stack_Ori=Image_Stack;
        Image_Stack_Ori=Image_Stack_Ori(:,:,size(Image_Stack_Ori,3):-1:1);
        if If_Hilbert == 0
            Image_Stack=2*real(ifft(Image_Stack_FFT.*repmat(Filter_Function,[size(Image_Stack_FFT,1) size(Image_Stack_FFT,2) 1]),[],3));
            Phase_Stack=angle(ifft(Image_Stack_FFT.*repmat(Filter_Function,[size(Image_Stack_FFT,1) size(Image_Stack_FFT,2) 1]),[],3));
            Phase_Stack=Phase_Stack(:,:,size(Phase_Stack,3):-1:1);

        else
            Image_Stack=2*abs(ifft(Image_Stack_FFT.*repmat(Filter_Function,[size(Image_Stack_FFT,1) size(Image_Stack_FFT,2) 1]),[],3));
            Phase_Stack=angle(ifft(Image_Stack_FFT.*repmat(Filter_Function,[size(Image_Stack_FFT,1) size(Image_Stack_FFT,2) 1]),[],3));
            Phase_Stack=Phase_Stack(:,:,size(Phase_Stack,3):-1:1);
        end
        Image_Stack=Image_Stack(:,:,size(Image_Stack,3):-1:1);

        %%
                %imagesc(squeeze(Image_Stack(:,6,:))')

    end
    %% Image Sharpness Estimation
    if If_Sharpness_Est == 1
        
        %% Post-FFT Axial Ave
        Reduced_Image=squeeze(mean(Image_Stack,2))';
       
        Temp=0;
        Length_Used_Post=floor(size(Reduced_Image,1)/AVE_postHilbert)*AVE_postHilbert;
        for p=1:AVE_postHilbert
           Temp=Temp+Reduced_Image((AVE_postHilbert-(p-1)):AVE_postHilbert:(Length_Used_Post-(p-1)),:);
        end
        Reduced_Image=Temp/AVE_postHilbert;   
        %% Generating Axial Profile for Reduced Image
        Z_Profile=Axial_Resolution*([0:(size(Reduced_Image,1)-1)]);
       Index_Range=[find(Z_Profile>Axial_ROI(1),1,'first') find(Z_Profile>Axial_ROI(2),1,'first')];
        
        
        %% Basic Index        
        [Gx, Gz]=gradient(Reduced_Image);
        if If_Lateral_Sharpness_Only ==1
            Gz=0;
        end
        S=sqrt(Gx.^2+Gz.^2)/mean(mean(Ave_Array))*4096;  %!!!!! �γ̭�lintensity norm, ���I��NN
        sharpness=mean(mean(S(Index_Range(1):Index_Range(2),:)));
        %% Intensity Profile
        intensity=mean(mean(Reduced_Image));
        I_Profile=mean(Reduced_Image,2);
        %% No normalization at all        
        S_NoNorm=S*mean(mean(Ave_Array))/4096;
        sharpness_NoNorm=mean(mean(S_NoNorm(Index_Range(1):Index_Range(2),:)));
        %% Depth Enhacement

        Linear_Enhace_Coef=[BI:(BF-BI)/(size(Reduced_Image,1)-1):BF]';
        Enhanced_Image=Reduced_Image.*repmat(Linear_Enhace_Coef,[1 size(Reduced_Image,2)]);
        Enhanced_Axial_Profile=mean(Enhanced_Image,2);
        %
        %
        %% Auto-brightness
        Auto_Brightness_TH=0.99;
        nbin=250;
        Histogram=hist(Enhanced_Image(:),nbin);
        Total=sum(Histogram);
        Int_Hist = cumsum(Histogram)/Total;

        C_max=find(Int_Hist>Auto_Brightness_TH,1,'first')/nbin*max(Enhanced_Image(:));
        C_min=C_max*0.1;
        Normalized_Image=(Enhanced_Image-C_min)/(C_max-C_min);
        Normalized_Image(Normalized_Image<0)=0;
        Normalized_Image(Normalized_Image>1)=1;
        %
        subplot(1,1,1)
        imagesc(Normalized_Image);
        colormap(gray);
        %% Processed Index        
        [Gx_Processed, Gz_Processed]=gradient(Normalized_Image);
        if If_Lateral_Sharpness_Only ==1
            Gz_Processed=0;
        end
        S_Processed=sqrt(Gx_Processed.^2+Gz_Processed.^2);
        sharpness_Processed=mean(mean(S_Processed(Index_Range(1):Index_Range(2),:)));
   
        %% Normalized Gradient
        [Gx_Norm_temp, Gz_Norm_temp]=gradient(Reduced_Image);
        if If_Lateral_Sharpness_Only ==1
            Gz_Norm_temp=0;
        end
        Gx_Norm=Gx_Norm_temp./Reduced_Image;
        Gz_Norm=Gz_Norm_temp./Reduced_Image;
        Gx_Norm(Reduced_Image==0)=0;
        Gz_Norm(Reduced_Image==0)=0;
        
        S_Norm=sqrt(Gx_Norm.^2+Gz_Norm.^2);    %�b��normalized S�U, noise���^�m�Q����
        %imagesc(abs(S_Norm));
        
        sharpness_Norm=mean(mean((S_Norm(Index_Range(1):Index_Range(2),:))));
        S_Profile=mean(S,2);
        S_Processed_Profile=mean(S_Processed,2);
        S_Norm_Profile=mean(S_Norm,2);
        
        subplot(4,1,1)
        plot(Z_Profile,S_Profile,'LineWidth',2)
        xlim([0 max(Z_Profile)])
        ylim([0 10])
        xlabel('Depth (micron)')
        ylabel('Sharpness')
        title(sprintf('Sharpness Profile, S=%.04f',sharpness))

        subplot(4,1,2)
        plot(Z_Profile,I_Profile,'LineWidth',2)
        xlim([0 max(Z_Profile)])
        ylim([0 30])
        xlabel('Depth (micron)')
        ylabel('Intensity')
        title(sprintf('Intensity Profile, I=%.04f',intensity))

        subplot(4,1,3)
        plot(Z_Profile,S_Processed_Profile,'LineWidth',2)
        xlim([0 max(Z_Profile)])
        ylim([0 0.15])
        xlabel('Depth (micron)')
        ylabel('Processed Sharpness')
        title(sprintf('Processed Sharpness Profile, S_p_r_o_c=%.04f',sharpness_Processed))
        
        subplot(4,1,4)
        plot(Z_Profile,S_Norm_Profile,'LineWidth',2)
        xlim([0 max(Z_Profile)])
        ylim([0 0.3])
        xlabel('Depth (micron)')
        ylabel('Normalized Sharpness')
        title(sprintf('Normalized Sharpness Profile, S_n_o_r_m=%.04f',sharpness_Norm))
        saveas(gcf,[Data_Save_Folder sprintf('%s_Sharpness.png',folder_list(QQQ).name)]);

        fprintf('\nSharpness Index: %g',sharpness);
        fprintf('\nProcessed Sharpness Index: %g',sharpness_Processed);
        fprintf('\nPixel-wise Normalized Sharpness Index: %g\n\n',sharpness_Norm);

        imwrite(Normalized_Image,[Data_Save_Folder sprintf('%s_Sharpness=%04g_Processed Sharpness=%04g_Pixel-wise Normalized Sharpness=%04g.png',folder_list(QQQ).name,sharpness,sharpness_Processed,sharpness_Norm)],'png');
        dlmwrite([Data_Save_Folder 'Sharpness.txt'],[QQQ sharpness sharpness_Processed sharpness_Norm sharpness_NoNorm intensity] ,'delimiter','\t','newline','pc','precision', '%.6f','-append');

    end
  
       %% Spectrogram for STFT
    if (If_STFT_1D == 1)
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
        Predicted_Central_Bin=nBin/2/Half_Axial_Spatial_Sampling_Frequency*Predicted_Central_Frequency;
        Max_Spectral_Frequency=3E8*Half_Axial_Spatial_Sampling_Frequency/2*1E6;   %/2 for roundtrip; *1E6 for micron to mm, TiSa ~=3.8E14
        Predicted_Spectral_Frequency=3E8/(Predicted_Central_Wavelength)*1E6;
        
        Spectral_Frequency_Array=0:(Max_Spectral_Frequency/(nBin/2)):Max_Spectral_Frequency;
        Spectral_Frequency_of_Interst_Array=3E8./Wavelength_of_Interest_Array*1E6;
        Ave_Spectrogram_Based_on_Wavelength=interp2(repmat(1:size(Ave_Spectrogram,2),[size(Ave_Spectrogram,1) 1]),repmat(Spectral_Frequency_Array',[1 size(Ave_Spectrogram,2)]),Ave_Spectrogram,repmat(1:size(Ave_Spectrogram,2),[length(Spectral_Frequency_of_Interst_Array) 1]),repmat(Spectral_Frequency_of_Interst_Array',[1 size(Ave_Spectrogram,2)]));
        [Max_Value Max_Index]=max(Ave_Spectrogram_Based_on_Wavelength,[],1);
        %%  17/11/12 More Data Analysis
        Spectrograqm_Contrast_Index=max(Ave_Spectrogram_Based_on_Wavelength,[],1)./min(Ave_Spectrogram_Based_on_Wavelength,[],1);
        plot(Spectrograqm_Contrast_Index);
        %%
        Ave_Spectrogram_Based_on_Wavelength_Norm=Ave_Spectrogram_Based_on_Wavelength./repmat(Max_Value,[size(Ave_Spectrogram_Based_on_Wavelength,1) 1]);
        
        Bin_FPC=[(1:size(Ave_Spectrogram_Based_on_Wavelength_Norm,2))' (Wavelength_of_Interest_Array(Max_Index)/Predicted_Central_Wavelength*N*ave_factor)'];
        
        dlmwrite([Data_Save_Folder sprintf('%s_Bin_FPC.txt',folder_list(QQQ).name)],Bin_FPC,'delimiter','\t','newline','pc','precision', '%.6f');
        subplot(1,1,1)
        imagesc(Ave_Spectrogram_Based_on_Wavelength_Norm);
        colormap(parula); %parula
        hold on
        plot(1:length(Max_Index),Max_Index,'ro');
        hold off
        set(gca,'YTick',(1:8:length(Wavelength_of_Interest_Array)),'YTickLabel',round(downsample(Wavelength_of_Interest_Array/Predicted_Central_Wavelength*N*ave_factor,8)*10)/10,'FontSize', 14);
        xlabel('Temporal (Axial) Bins');
        ylabel('Frame per Carrier');
        pbaspect([2 1 1]);
        saveas(gcf,[Data_Save_Folder sprintf('%s_FPC.png',folder_list(QQQ).name)]);
        
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

        BI=-0.1;
        BF=1;
        Linear_Enhace_Coef=reshape([BI:(BF-BI)/(size(Image_Stack,3)-1):BF],1,1,[]);

            Image_Stack=Image_Stack.*repmat(Linear_Enhace_Coef,[size(Image_Stack,1) size(Image_Stack,2) 1]);

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
        
        %%
        if If_Partial_Spectrum == 0
            imwrite(Npoint_Image,[Data_Save_Folder sprintf('%s_Npoint.png',folder_list(QQQ).name)],'png');
            imwrite(uint16(Npoint_Image*2^16),[Data_Save_Folder sprintf('%s_Npoint.tif',folder_list(QQQ).name)],'tiff');

        else
            imwrite(Npoint_Image,[Data_Save_Folder sprintf('%s_Npoint_Filtered_BW%g_Offset%g.png',folder_list(QQQ).name,BW,Offset)],'png');
            imwrite(uint16(Npoint_Image*2^16),[Data_Save_Folder sprintf('%s_Npoint_Filtered_BW%g_Offset%g.tif',folder_list(QQQ).name,BW,Offset)],'tiff');
        end

    end
    
    %% Post N-point Speckle Axial
    
    %% Save Raw File
    if If_Save_Raw == 1
        if If_Decorrelation == 0
            Processed_Data_Path=[Data_Save_Folder '\' folder_list(QQQ).name '\' sprintf('%s.raw',folder_list(QQQ).name)];
            mkdir([Data_Save_Folder '\' folder_list(QQQ).name '\']);
        else
            Processed_Data_Path=[Data_Save_Folder '\' folder_list(QQQ).name '\' sprintf('%s_Decorr_Pre%g_Post%g.raw',folder_list(QQQ).name,AVE_preDecorr,AVE_postDecorr)];            
        end
        Image_Stack=Image_Stack(:,:,size(Image_Stack,3):-1:1);
        fid = fopen(Processed_Data_Path, 'w+');
        fwrite(fid, Image_Stack, 'double');
        %fwrite(fid, Phase_Stack, 'double');
        
        fclose(fid);
        fclose('all');
    end
end



