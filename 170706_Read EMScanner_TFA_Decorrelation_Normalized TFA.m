clear all
%% Operations
If_DF=1;
If_FF=0;
If_normalization=1;
Normalization_ROI=[400 600];

If_Partial_Spectrum=0;
BW=0.005;        %Guassian sig
Offset=0;    %percent

If_Decorrelation=0; %Along Z, ~= along Time
AVE_preDecorr=2;
AVE_postDecorr=16/AVE_preDecorr;
Decorr_Skip=8;


If_Npoint=0; %Along Z, ~= along Time
N=4;
AVE_preNp=2;
AVE_postNp=16/AVE_preDecorr/N;
Auto_Brightness_TH=0.99;

If_Decorrelation_Y=0;
AVE_postDecorr_Y=1;

If_Save_Raw=0;

If_STFT_1D=1;   % only use the central 1D map for STFT, reducing data szie
Time_Window=500;
nBin=1024;
Spectrogram_X_ROI=[200 800];
Spectrogram_Y_ROI=[1 11];
Wavelength_of_Interest_Array=0.6:0.01:1.2; %micron
RI=1.406;    %Refractive Index
Predicted_Central_Wavelength=0.78;
%


root_folder_path='C:\TuanShu\170706_PZT Scanning Range Test\';
last_folder_name='TiO2_sequential8';
Number=[];

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


If_AlterAxialProcess=1;

% Data format related
Row=1024;
Colomn=11;
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
    folder_path=[parent_folder_path '\' folder_list(QQQ).name '\'];
    cd(folder_path);
    micron_per_frame=0.2/ave_factor/N;
    Offset_1=0;
    Offset_2=0; 

    Max_Number_of_Frame=3000000;

    %%

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
    %%
    file_list=downsample(file_list,floor(Product_of_Axial_Decimation_Factor_and_Ave_Factor/ave_factor),Phase_Shift);
    Frame_Number_in_File_List=downsample(Frame_Number_in_File_List,floor(Product_of_Axial_Decimation_Factor_and_Ave_Factor/ave_factor),Phase_Shift);
    Frame=length(file_list);
    
    %%
    Image_Stack=zeros(Row,Colomn,Frame);
    for p=1:length(file_list)

        %file_path=[folder_path file_list(p).name];

        fin=fopen(file_list(p).name);

        fseek(fin, Byte_Skip+Row*Colomn*2*(Frame_Number_in_File_List(p)-1), 'bof');
        %QQ=fread(fin,[Row,Colomn],'uint16');   
        Image_Stack(:,:,p)=fread(fin,[Row,Colomn],'uint16'); %*Frame   不知為何, 看起來就像是要除16
        %%
        %fclose(fin);
        fclose('all');
        disp(p);
    end
        
        %error('Multiple Files in Folder');
    
    
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
    
    if If_Partial_Spectrum ==1
        Image_Stack_FFT=fft(Image_Stack,[],3);
        %%
        LB=0.01;         % ratio
        RB=0.5;
        Image_Stack_FFT(:,:,1:round(size(Image_Stack_FFT,3)*LB))=0;
        %Image_Stack_FFT(:,:,round(size(Image_Stack_FFT,3)/2):end)=0;
        Image_Stack_FFT(:,:,round(size(Image_Stack_FFT,3)*RB):end)=0;
        
        Mean_Spectrum=squeeze(mean(mean(abs(Image_Stack_FFT),1),2));    %%only used for center finding
        Mean_Spectrum(Mean_Spectrum<max(Mean_Spectrum)*0.2)=0;          %only >20% peak are considered
        Peak_Index=sum([1:length(Mean_Spectrum)]'.*Mean_Spectrum)/sum(Mean_Spectrum);
        Filter_Function=zeros(1,1,length(Mean_Spectrum));
        Filter_Function(1,1,:)=gaussmf([1:length(Mean_Spectrum)]',[BW*length(Mean_Spectrum) Peak_Index+round(size(Image_Stack_FFT,3)*Offset)]);
        plot(squeeze(mean(mean(abs(Image_Stack_FFT),1),2)).*squeeze(Filter_Function));
        %plot(squeeze(mean(mean(abs(Image_Stack_FFT),1),2)))
%%
        %NNN=101;
        %plot(squeeze(abs(Image_Stack_FFT(NNN,6,:))));
        %%
        Image_Stack=2*real(ifft(Image_Stack_FFT.*repmat(Filter_Function,[size(Image_Stack_FFT,1) size(Image_Stack_FFT,2) 1]),[],3));
        %Image_Stack=2*abs(ifft(Image_Stack_FFT,[],3));
        %Image_Stack=2*real(ifft(Image_Stack_FFT,[],3));
    end
    
    
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
        imwrite(Decorr_Image,[Data_Save_Folder sprintf('%s_Decorr_Pre%g_Post%g_Skip%g.png',folder_list(QQQ).name,AVE_preDecorr,AVE_postDecorr,Decorr_Skip)],'png');
    end
    
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
        
        Npoint_Image=squeeze(mean(Image_Stack,2))';
        Npoint_Image=Npoint_Image(size(Npoint_Image,1):-1:1,:);
        
        nbin=250;
        Histogram=hist(Npoint_Image(:),nbin);
        Total=sum(Histogram);
        Int_Hist = cumsum(Histogram)/Total;

        C_max=find(Int_Hist>Auto_Brightness_TH,1,'first')/nbin*max(Npoint_Image(:));
        C_min=C_max*0.05;

        Npoint_Image=(Npoint_Image-C_min)/(C_max-C_min);
        Npoint_Image(Npoint_Image>1)=1;
        Npoint_Image(Npoint_Image<0)=0;
        if If_Partial_Spectrum == 0
            imwrite(Npoint_Image,[Data_Save_Folder sprintf('%s_Npoint.png',folder_list(QQQ).name)],'png');
            imwrite(uint16(Npoint_Image*2^16),[Data_Save_Folder sprintf('%s_Npoint.tif',folder_list(QQQ).name)],'tiff');

        else
            imwrite(Npoint_Image,[Data_Save_Folder sprintf('%s_Npoint_Filtered_BW%g_Offset%g.png',folder_list(QQQ).name,BW,Offset)],'png');
            imwrite(uint16(Npoint_Image*2^16),[Data_Save_Folder sprintf('%s_Npoint_Filtered_BW%g_Offset%g.tif',folder_list(QQQ).name,BW,Offset)],'tiff');
        end
        imagesc(Npoint_Image);
        colormap(gray);
    end
    
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
        imwrite(Decorr_Image,[Data_Save_Folder sprintf('%s_Decorr_Y.png',folder_list(QQQ).name)],'png');
    end
    
    if If_STFT_1D
        Image_Virtual_B_Scan=reshape(Image_Stack(Spectrogram_X_ROI(1):Spectrogram_X_ROI(2),Spectrogram_Y_ROI(1):Spectrogram_Y_ROI(2),:),(Spectrogram_X_ROI(2)-Spectrogram_X_ROI(1)+1)*(Spectrogram_Y_ROI(2)-Spectrogram_Y_ROI(1)+1),[]);
        Image_Spectrogram=zeros(size(Image_Virtual_B_Scan,1),nBin/2+1,floor(size(Image_Virtual_B_Scan,2)/Time_Window));
        for p=1:size(Image_Virtual_B_Scan,1)
            Image_Spectrogram(p,:,:)=spectrogram(Image_Virtual_B_Scan(p,:),Time_Window,0,nBin);  %1st 0 for zero overlapping, nBin for DFT bins, it will outout only half spectrum (512+1)
            disp(p);
        end
        Ave_Spectrogram=squeeze(mean(abs(Image_Spectrogram(Spectrogram_X_ROI(1):Spectrogram_X_ROI(2),:,:)),1));
        imagesc(abs(Ave_Spectrogram));  %max frequency = 0.5 * sampling frequency

        Axial_Spatial_Sampling_Resolution=Predicted_Central_Wavelength/2/RI/N/ave_factor;                    %micron
        Axial_Spatial_Sampling_Resolution_OPD=Predicted_Central_Wavelength/2/N/ave_factor;                    %micron
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
    if If_Save_Raw == 1
        if If_Decorrelation == 0
            Processed_Data_Path=[Data_Save_Folder sprintf('%s.raw',folder_list(QQQ).name)];
        else
            Processed_Data_Path=[Data_Save_Folder sprintf('%s_Decorr_Pre%g_Post%g.raw',folder_list(QQQ).name,AVE_preDecorr,AVE_postDecorr)];            
        end

        fid = fopen(Processed_Data_Path, 'w+');
        fwrite(fid, flip(Image_Stack,3), 'double');
        fclose(fid);
        fclose('all');
    end
end
%%
% NNN=656;
% imagesc(squeeze(abs(Image_Spectrogram(NNN,:,:))));  %lateral is time
% 

%%

%Ave_Spectrogram_Based_on_Wavelength=interp2(repmat(Spectral_Frequency_Array',[1 size(Ave_Spectrogram,2)]),repmat(1:size(Ave_Spectrogram,2),[size(Ave_Spectrogram,1) 1]),Ave_Spectrogram,repmat(Spectral_Frequency_of_Interst_Array',[1 size(Ave_Spectrogram,2)]),repmat(1:size(Ave_Spectrogram,2),[size(Ave_Spectrogram,1) 1]));
%%


