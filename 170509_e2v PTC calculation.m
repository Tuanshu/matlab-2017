clear all
%%
root_folder_path='D:\170509_OctoPlus_PTC test\Diff Current\';
Data_Name='170509_OctoPlus_PTC test_DiffCurrent';
Current_Array=[150 200 250 300 400 500 600 700 800 900 1000];


Horizontal_Size=2048;
Size=25;

H_Offset=501;
H_Range=1000;
Horizontal_ROI=[H_Offset H_Offset+H_Range];

Noise_Removal_Method=2; %1: do nothing, 2: normalization, 3: filtering, 4: 4-point

% Feature for Filtering mode
Filtering_Ratio=0.02;       %意思是, 若Array長度為A, 則將外側(低頻)之A*Filtering_Ratio個pixels設為0
% Feature for 4-point mode
Averaging_Factor=1;
N=4;          %for N-point
Unbias_Estimator=0.9213;

Mean_Array_Map=zeros(Horizontal_ROI(2)-Horizontal_ROI(1)+1,length(Current_Array));
STD_Array_Map=zeros(Horizontal_ROI(2)-Horizontal_ROI(1)+1,length(Current_Array));

Mean_Array_for_norm=zeros(1,length(Current_Array));

for p=1:length(Current_Array)
    file_name=sprintf('%g',Current_Array(p));
    folder_path=[root_folder_path '\'];
    file_path=[folder_path file_name];

    Data_Save_Folder='D:\Processed Data\';

    Processed_Data_Path=[Data_Save_Folder Data_Name '_' file_name '.bin'];

    cd(folder_path);

    Data_Image=double(imread(file_path,'tiff'));

    Number_of_Filter_Pixel_one_side=size(Data_Image,1)*Filtering_Ratio/2;   %總共有左右兩個side
    
    Mean_Array_for_norm=mean(Data_Image,2); % 固定用全張影像做, 有好有壞, 優點: 不會因為sample點數太少而不準, 缺點: 若影像有部分飽和則不準
    All_Mean=mean(Mean_Array_for_norm);
    Mean_Image_for_norm=repmat(Mean_Array_for_norm,[1 size(Data_Image,2)]);
    Data_Image_norm=Data_Image./Mean_Image_for_norm*All_Mean;
    
    Data_Image_FFT=fft((Data_Image-repmat(mean(Data_Image,1),[size(Data_Image,1) 1])),[],1);
    
    Data_Image_FFT_Filter=Data_Image_FFT;
    Data_Image_FFT_Filter(1:Number_of_Filter_Pixel_one_side,:)=0;
    Data_Image_FFT_Filter((size(Data_Image_FFT_Filter,1)-Number_of_Filter_Pixel_one_side+1):size(Data_Image_FFT_Filter,1),:)=0;
    Data_Image_Filter=real(ifft(Data_Image_FFT_Filter,[],1))+repmat(mean(Data_Image,1),[size(Data_Image,1) 1]);

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
    
    
    if Noise_Removal_Method == 1
        Mean_Array_Map(Horizontal_ROI(1):Horizontal_ROI(2),p)=mean(Data_Image(:,Horizontal_ROI(1):Horizontal_ROI(2)),1);
        STD_Array_Map(Horizontal_ROI(1):Horizontal_ROI(2),p)=std(Data_Image(:,Horizontal_ROI(1):Horizontal_ROI(2)),1,1);
    elseif Noise_Removal_Method == 2
        Mean_Array_Map(Horizontal_ROI(1):Horizontal_ROI(2),p)=mean(Data_Image_norm(:,Horizontal_ROI(1):Horizontal_ROI(2)),1);
        STD_Array_Map(Horizontal_ROI(1):Horizontal_ROI(2),p)=std(Data_Image_norm(:,Horizontal_ROI(1):Horizontal_ROI(2)),1,1);    
    elseif Noise_Removal_Method == 3
        Mean_Array_Map(Horizontal_ROI(1):Horizontal_ROI(2),p)=mean(Data_Image_Filter(:,Horizontal_ROI(1):Horizontal_ROI(2)),1);
        STD_Array_Map(Horizontal_ROI(1):Horizontal_ROI(2),p)=std(Data_Image_Filter(:,Horizontal_ROI(1):Horizontal_ROI(2)),1,1);        
    elseif Noise_Removal_Method == 4
        Mean_Array_Map(Horizontal_ROI(1):Horizontal_ROI(2),p)=mean(Data_Image(:,Horizontal_ROI(1):Horizontal_ROI(2)),1);
        STD_Array_Map(Horizontal_ROI(1):Horizontal_ROI(2),p)=mean(Data_Image_Npoint(:,Horizontal_ROI(1):Horizontal_ROI(2)),1)*Unbias_Estimator;        
    end
    disp(p);
end
All_Meam=mean(Mean_Array_for_norm);

% %%
% NN=6;
% 
% plot(real(Data_Image_FFT(:,NN)));
% 

%


Mean_Array_Map_1D=Mean_Array_Map(:);
STD_Array_Map_1D=STD_Array_Map(:);
VAR_Array_Map_1D=STD_Array_Map_1D.^2;
scatter(Mean_Array_Map_1D,VAR_Array_Map_1D,Size,'filled');
xlabel('Signal (DN)','fontsize',15);
ylabel('Variance (DN^2)','fontsize',15);
set(gca,'fontsize',15)
ylim([0 130]);
xlim([0 4096]);
saveas(gcf,[Data_Save_Folder Data_Name sprintf('_Method_%g',Noise_Removal_Method) '.png']);
