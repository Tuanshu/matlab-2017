clear all
%%
disc='D:\';
root_folder_path='170414_OctoPlus+500proj';
mid_folder_path='';
sub_folder_path='Forearm';
Number=3;


C_max=72;%48;   %32
C_min=0;%0.8;  %0.8

%root_folder_path='I:\161223_e2v EM1\CeYAG\Different Exp Time';
Averaging_Factor=32;
Product_of_Axial_Decimation_Factor_and_Ave_Factor=32;
Vertical_Maximum_Depth=10000;    %pixel 1285

folder_path=[disc root_folder_path '\' mid_folder_path '\' sub_folder_path '_' sprintf('%g',Number) '\'];

Resample_Ratio=1;

Data_Name=sub_folder_path;
it=29;
Gain=-24;
Display_Gain=4;
Horizontal_Size=2048;
Size=25;


SR_Lateral=0.18;
SR_Axial=0.279;
Row_Binning_Factor=5;
Axial_ave_Factor=1;
%

H_Offset_View=1; %80  200
H_Range_View=2000; %750
Horizontal_ROI_View=[H_Offset_View H_Offset_View+H_Range_View];

%
H_Offset=480;
H_Range=250;
Horizontal_ROI=[H_Offset H_Offset+H_Range];

If_norm=1; %0: do nothing, 1: normalization

% Feature for Filtering mode
Filtering_Ratio=0.02;       %意思是, 若Array長度為A, 則將外側(低頻)之A*Filtering_Ratio個pixels設為0
% Feature for 4-point mode
N=4;          %for N-point
Unbias_Estimator=0.9213;
% 
% Mean_Array_Map=zeros(Horizontal_ROI(2)-Horizontal_ROI(1)+1,length(Current_Array));
% STD_Array_Map=zeros(Horizontal_ROI(2)-Horizontal_ROI(1)+1,length(Current_Array));
% 
% Mean_Array_for_norm=zeros(1,length(Current_Array));

Data_Save_Folder=[disc root_folder_path '\' mid_folder_path '\Processed Data\'];%'C:\Users\Owner\Desktop\031317\Processed Data\';
mkdir(Data_Save_Folder);


Current_Array=[1];%[3:3:24];

file_list=dir(folder_path);
file_list=file_list(4:end);
Stitched_Image=[];
for p=1:length(file_list)
    file_name=file_list(p).name;
    file_path=[folder_path '\' file_name];


    %Processed_Data_Path=[Data_Save_Folder Data_Name '_' file_name '.bin'];

    %cd(folder_path);

    Data_Image_Raw=double(imread(file_path,'tiff'));
    Data_Image_Raw_Decimated=Data_Image_Raw(1:floor(Product_of_Axial_Decimation_Factor_and_Ave_Factor/Averaging_Factor):end);
       
    if If_norm == 0
           Data_Image=Data_Image_Raw;
    elseif If_norm == 1
        Mean_Array_for_norm=mean(Data_Image_Raw(:,Horizontal_ROI(1):Horizontal_ROI(2)),2); % 改成可調ROI
        All_Mean=mean(Mean_Array_for_norm);
        Mean_Image_for_norm=repmat(Mean_Array_for_norm,[1 size(Data_Image_Raw,2)]);
        Data_Image=Data_Image_Raw./Mean_Image_for_norm*All_Mean;
    end
    
    Averaged_Length=floor(size(Data_Image,1)/Averaging_Factor);
    Temp=0;
    for q=1:Averaging_Factor
       Temp=Temp+Data_Image((Averaging_Factor-(q-1)):Averaging_Factor:(Averaging_Factor*Averaged_Length)-(q-1),:);
    end
    Data_Image_Ave=Temp/Averaging_Factor;
    Data_Image_Npoint=zeros(floor(size(Data_Image_Ave,1)/N),size(Data_Image_Ave,2));
    for q=1:size(Data_Image_Npoint,1)
        Temp=Data_Image_Ave(((q-1)*N+1):(q*N),:);
        Data_Image_Npoint(q,:)=((N*sum(Temp.^2,1)-sum(Temp,1).^2).^0.5)*(2^0.5)/N;
    end
    Stitched_Image=[Stitched_Image; Data_Image_Npoint];
end
Stitched_Image=Stitched_Image(size(Stitched_Image,1):-1:1,:);
%imagesc(Data_Image_Npoint);
%
if size(Stitched_Image,1)>Vertical_Maximum_Depth
    Stitched_Image=Stitched_Image(1:Vertical_Maximum_Depth,:);
end

% SNR calc
Safe_Margin=10;
Blur=5;
Blur_Filter=fspecial('gaussian',Blur,Blur/3);
Axial_Profile_Withglass=mean(Stitched_Image(:,Horizontal_ROI_View(1):Horizontal_ROI_View(2)),2);
%plot(Axial_Profile);
[Glass_Interface_value Glass_Interface_index]=max(Axial_Profile_Withglass);

Data_Image_Noglass=Stitched_Image((Glass_Interface_index+Safe_Margin):end,Horizontal_ROI_View(1):Horizontal_ROI_View(2));
Data_Image_Beforeglass=Stitched_Image(1:(Glass_Interface_index-Safe_Margin),Horizontal_ROI_View(1):Horizontal_ROI_View(2));

Noise=std(Data_Image_Beforeglass(:));

Data_Image_Noglass_Blur=filter2(Blur_Filter,Data_Image_Noglass);
Axial_Profile_Noglass=mean(Data_Image_Noglass,2);

Signal=max(Data_Image_Noglass_Blur(:));

Histogram_Withglass=hist(Stitched_Image(:),[0:0.01:200]);
Histogram_Noglass=hist(Data_Image_Noglass(:),[0:0.01:200]);

%plot([0:0.01:200],Histogram_Withglass,[0:0.01:200],Histogram_Noglass);

SNR=20*log10(Signal/Noise);

fprintf('\nNoise=%.3f, Signal=%.3f, SNR=%.3f dB\n',Noise,Signal,SNR)

%


%
Stitched_Image_FOV=Stitched_Image(:,Horizontal_ROI_View(1):Horizontal_ROI_View(2));
Data_Image_Npoint_normalized_BeforeLB=Display_Gain*(Stitched_Image_FOV-C_min)/(C_max-C_min);

    % Row Binning
Temp=0;
Lateral_Length_Original=size(Data_Image_Npoint_normalized_BeforeLB,2);
Reduced_Length=Lateral_Length_Original/Row_Binning_Factor;
for p=1:Row_Binning_Factor
   Temp=Temp+Data_Image_Npoint_normalized_BeforeLB(:,(Row_Binning_Factor-(p-1)):Row_Binning_Factor:(Row_Binning_Factor*Reduced_Length)-(p-1));
end
Data_Image_Npoint_normalized=Temp/Row_Binning_Factor;




Temp=0;
Axial_Length_Original=size(Data_Image_Npoint_normalized,1);
Axial_Length_Used=floor(Axial_Length_Original/Axial_ave_Factor)*Axial_ave_Factor;
Reduced_Length=Axial_Length_Used/Axial_ave_Factor;
for p=1:Axial_ave_Factor
   Temp=Temp+Data_Image_Npoint_normalized((Axial_ave_Factor-(p-1)):Axial_ave_Factor:(Axial_ave_Factor*Reduced_Length)-(p-1),:);
end
Data_Image_Reduced=Temp/Axial_ave_Factor;


Data_Image_Npoint_normalized_Resampled=zeros(size(Data_Image_Npoint_normalized,1),size(Data_Image_Npoint_normalized,2)*Resample_Ratio);
for p=1:Resample_Ratio
    Data_Image_Npoint_normalized_Resampled(:,p:Resample_Ratio:Resample_Ratio*(size(Data_Image_Npoint_normalized,2)-1)+p)=Data_Image_Npoint_normalized;
end

Data_Image_Npoint_normalized(Data_Image_Npoint_normalized<0)=0;
Data_Image_Npoint_normalized(Data_Image_Npoint_normalized>1)=1;

Data_Image_Reduced(Data_Image_Reduced<0)=0;
Data_Image_Reduced(Data_Image_Reduced>1)=1;


Data_Image_Npoint_normalized_Resampled(Data_Image_Npoint_normalized_Resampled<0)=0;
Data_Image_Npoint_normalized_Resampled(Data_Image_Npoint_normalized_Resampled>1)=1;




subplot(1,1,1)
imagesc(Data_Image_Npoint_normalized_Resampled);
%axis off
axis equal
caxis([0 1]);
X_Label_Array=[0  100 200 300 400 500 600];
X_Label_Index_Array=X_Label_Array/(SR_Lateral)*Resample_Ratio;
xlim([0 size(Data_Image_Npoint_normalized_Resampled,2)]);
Y_Label_Array=[0 50 100 150 200 250 300 350 400];
Y_Label_Index_Array=round(Y_Label_Array/SR_Axial);
set(gca,'xtick',X_Label_Index_Array,'xticklabel',X_Label_Array,'ytick',Y_Label_Index_Array,'yticklabel',Y_Label_Array);
ylim([0 size(Data_Image_Npoint_normalized_Resampled,1)]);
xlabel('Lateral Position (micron)');
ylabel('Axial Position (micron)');

colormap(gray);

if If_norm ==0
    imwrite(Data_Image_Npoint_normalized,[Data_Save_Folder sprintf('%s_%s_%d_N%.1f_S%.1f_SNR%.1fdB',root_folder_path,sub_folder_path,Number,Noise,Signal,SNR) '.png'],'png');
    saveas(gcf,[Data_Save_Folder sprintf('%s_%s_%d_N%.2f_S%.2f_SNR%.2fdB',root_folder_path,sub_folder_path,Number,Noise,Signal,SNR) '_ScaleBar.png']);
    imwrite(Data_Image_Npoint_normalized_Resampled,[Data_Save_Folder sprintf('%s_%s_%d_N%.2f_S%.2f_SNR%.2fdB',root_folder_path,sub_folder_path,Number,Noise,Signal,SNR) '_Resampled.png'],'png');
    imwrite(Data_Image_Reduced,[Data_Save_Folder sprintf('%s_%s_%d_N%.2f_S%.2f_SNR%.2fdB',root_folder_path,sub_folder_path,Number,Noise,Signal,SNR) '_Reduced.png'],'png');

    fid = fopen([Data_Save_Folder sprintf('%s_%s_%d_N%.2f_S%.2f_SNR%.2fdB',root_folder_path,sub_folder_path,Number,Noise,Signal,SNR) '.raw'], 'w+');
    fwrite(fid, Stitched_Image_FOV', 'single');
    fclose(fid);
    
elseif If_norm ==1
    imwrite(Data_Image_Npoint_normalized,[Data_Save_Folder sprintf('%s_%s_%d_N%.2f_S%.2f_SNR%.2fdB',root_folder_path,sub_folder_path,Number,Noise,Signal,SNR) '_Norm.png'],'png');
    saveas(gcf,[Data_Save_Folder sprintf('%s_%s_%d_N%.2f_S%.2f_SNR%.2fdB',root_folder_path,sub_folder_path,Number,Noise,Signal,SNR) '_ScaleBar.png']);
    imwrite(Data_Image_Npoint_normalized_Resampled,[Data_Save_Folder sprintf('%s_%s_%d_N%.2f_S%.2f_SNR%.2fdB',root_folder_path,sub_folder_path,Number,Noise,Signal,SNR) '_Norm_Resampled.png'],'png');
    imwrite(Data_Image_Reduced,[Data_Save_Folder sprintf('%s_%s_%d_N%.2f_S%.2f_SNR%.2fdB',root_folder_path,sub_folder_path,Number,Noise,Signal,SNR) '_Norm_Reduced.png'],'png');

    fid = fopen([Data_Save_Folder sprintf('%s_%s_%d_N%.2f_S%.2f_SNR%.2fdB',root_folder_path,sub_folder_path,Number,Noise,Signal,SNR) '_Norm.raw'], 'w+');
    fwrite(fid, Stitched_Image_FOV', 'single');
    fclose(fid);
end
    
