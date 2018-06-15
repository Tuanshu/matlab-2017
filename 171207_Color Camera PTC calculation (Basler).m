clearvars 
%%
cd('D:\MATLAB\TsLib');
%%
Max_Number_of_Frame=50;
X_ROI=[1201 1300]-250;
Y_ROI=[951 1050]+400;
Number=[1];%9
Set='Basler';
Color='G';
%% Generating Bayer Mask
X_Offset=0;
Y_Offset=0;
Mask=zeros(X_ROI(2)-X_ROI(1)+1,Y_ROI(2)-Y_ROI(1)+1);
X_Grid=repmat([1:(X_ROI(2)-X_ROI(1)+1)]',[1 Y_ROI(2)-Y_ROI(1)+1])+X_Offset;
Y_Grid=repmat([1:(Y_ROI(2)-Y_ROI(1)+1)],[X_ROI(2)-X_ROI(1)+1 1])+X_Offset;
R_Mask=(mod(X_Grid,2)==0) & (mod(Y_Grid,2)==0);
B_Mask=(mod(X_Grid+1,2)==0) & (mod(Y_Grid+1,2)==0);
G_Mask=~R_Mask & ~B_Mask;
imagesc(G_Mask);

R_Mask_1D=R_Mask(:);
G_Mask_1D=G_Mask(:);
B_Mask_1D=B_Mask(:);

if strcmp(Color,'R')
    Mask_1D=R_Mask_1D(:);
elseif strcmp(Color,'G')
    Mask_1D=G_Mask_1D(:);
elseif strcmp(Color,'B')
    Mask_1D=B_Mask_1D(:);
end
Mask_1D=double(Mask_1D);
Mask_1D(Mask_1D==0)=NaN;
%% PTC related param
Noise_Removal_Method=1; %1: do nothing, 2: normalization, 3: filtering, 4: 4-point

% Feature for Filtering mode
Filtering_Ratio=0.02;       %意思是, 若Array長度為A, 則將外側(低頻)之A*Filtering_Ratio個pixels設為0
% Feature for 4-point mode
Averaging_Factor=1;
N=4;          %for N-point
Unbias_Estimator=0.9213;

if strcmp(Set,'Basler')
    root_folder_path='D:\';
    last_folder_name='171207_DiMirau Test';
    Row=2448;
    Colomn=2048;
    Bit_per_Pixel=24;
end
Frame_Size=Row*Colomn*Bit_per_Pixel/8;
   %Header=fread(fin,Header_Size,'ulong');

%fseek(fin, 256*4, 'bof');


% Scan Parent Folder
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
Matrix_Record=[];
Axial_Profile_Matrix=[];
Axial_Profile_10Perc_Matrix=[];
Noise_Array=[];
for QQQ=1:length(folder_list)

    folder_path=[parent_folder_path '\' folder_list(QQQ).name '\'];

    %cd(folder_path);

    %%
    file_list_ori=dir(folder_path);
    Original_Length_file_list=length(file_list_ori);
    if Original_Length_file_list <3
        continue;
    end
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
        Frame_Number_Calculated_Based_on_Size=(file_list_ori(p).bytes)/Frame_Size;
        if Frame_Number_Calculated_Based_on_Size>1.1
            file_list=[file_list;repmat(file_list_ori(p),[Frame_Number_Calculated_Based_on_Size 1])]; 
        else
            file_list=[file_list file_list_ori(p)];
        end
        Frame_Number_in_File_List=[Frame_Number_in_File_List;[1:Frame_Number_Calculated_Based_on_Size]']; 
    end

    %%
    Frame=min(length(file_list),Max_Number_of_Frame);
    %%
    X=[1:Frame];
    Image_Stack=zeros(X_ROI(2)-X_ROI(1)+1,Y_ROI(2)-Y_ROI(1)+1,Frame);
    for p=1:Frame
        file_path=[folder_path file_list(p).name];
        %QQ=fread(fin,[Row,Colomn],'uint16');
        if  Bit_per_Pixel == 8
            fin=fopen(file_path);
            fseek(fin, HeaderSize+Frame_Size*(Frame_Number_in_File_List(p)-1), 'bof');
            Image_Stack(:,:,p)=fread(fin,[Row,Colomn],'uint8');
        elseif Bit_per_Pixel == 16
            fin=fopen(file_path);
            fseek(fin, HeaderSize+Frame_Size*(Frame_Number_in_File_List(p)-1), 'bof');
            Image_Stack(:,:,p)=fread(fin,[Row,Colomn],'uint16');
        elseif Bit_per_Pixel == 24 %Color, RGB24
            Image_Temp=imread(file_path,'tiff');
            if strcmp(Color,'R')
                Image_Stack(:,:,p)=Image_Temp(Y_ROI(1):Y_ROI(2),X_ROI(1):X_ROI(2),1)';
            elseif strcmp(Color,'G')
                Image_Stack(:,:,p)=Image_Temp(Y_ROI(1):Y_ROI(2),X_ROI(1):X_ROI(2),2)';
            elseif strcmp(Color,'B')
                Image_Stack(:,:,p)=Image_Temp(Y_ROI(1):Y_ROI(2),X_ROI(1):X_ROI(2),3)';
            elseif strcmp(Color,'M')
                Image_Stack(:,:,p)=mean(Image_Temp(Y_ROI(1):Y_ROI(2),X_ROI(1):X_ROI(2),:),3)';
            end
        else
            continue;
        end
        fclose('all');
            %%
        disp(p);
    end
    %% PTC related Calcuation
    
    Number_of_Filter_Pixel_one_side=size(Image_Stack,3)*Filtering_Ratio/2;   %總共有左右兩個side

    Mean_Array_for_norm=mean(mean(Image_Stack,2),1); % 固定用全張影像做, 有好有壞, 優點: 不會因為sample點數太少而不準, 缺點: 若影像有部分飽和則不準
    Mean_Stack_for_norm=repmat(Mean_Array_for_norm,[size(Image_Stack,1) size(Image_Stack,2) 1]);
    
    Averaged_Length=floor(size(Image_Stack,3)/Averaging_Factor);
    Temp=0;
    for q=1:Averaging_Factor
       Temp=Temp+Image_Stack(:,:,(Averaging_Factor-(q-1)):Averaging_Factor:(Averaging_Factor*Averaged_Length)-(q-1));
    end
    Image_Stack_Ave=Temp/Averaging_Factor;
    if Noise_Removal_Method == 1
        Mean_Array_Map(:,:,QQQ)=mean(Image_Stack(:,:,:),3);
        STD_Array_Map(:,:,QQQ)=std(Image_Stack(:,:,:),0,3);
    elseif Noise_Removal_Method == 2
        All_Mean=mean(Mean_Array_for_norm);
        Image_Stack_norm=Image_Stack./Mean_Stack_for_norm*All_Mean;
        Mean_Array_Map(:,:,QQQ)=mean(Image_Stack_norm(:,:,:),3);
        STD_Array_Map(:,:,QQQ)=std(Image_Stack_norm(:,:,:),0,3);    
    elseif Noise_Removal_Method == 3
        Image_Stack_FFT=fft((Image_Stack-Mean_Stack_for_norm),[],3);
        Image_Stack_FFT_Filter=Image_Stack_FFT;
        Image_Stack_FFT_Filter(:,:,1:Number_of_Filter_Pixel_one_side)=0;
        Image_Stack_FFT_Filter(:,:,(size(Image_Stack_FFT_Filter,1)-Number_of_Filter_Pixel_one_side+1):size(Image_Stack_FFT_Filter,1))=0;
        Image_Stack_Filter=real(ifft(Image_Stack_FFT_Filter,[],3))+Mean_Stack_for_norm;

        Mean_Array_Map(:,:,QQQ)=mean(Image_Stack_Filter(:,:,:),3);
        STD_Array_Map(:,:,QQQ)=std(Image_Stack_Filter(:,:,:),0,3);        
    elseif Noise_Removal_Method == 4
        
        Image_Stack_Npoint=zeros(size(Image_Stack_Ave,1),size(Image_Stack_Ave,2),floor(size(Image_Stack_Ave,3)/N));
        for q=1:size(Image_Stack_Npoint,3)
            Temp=Image_Stack_Ave(:,:,((q-1)*N+1):(q*N));
            Image_Stack_Npoint(:,:,q)=((N*sum(Temp.^2,3)-sum(Temp,3).^2).^0.5)*(2^0.5)/N;
        end
        Mean_Array_Map(:,:,QQQ)=mean(Image_Stack(:,:,:),3);
        STD_Array_Map(:,:,QQQ)=mean(Image_Stack_Npoint(:,:,3),3)*Unbias_Estimator;        
    end
    disp(p);
    
end
 

%%
fclose('all');

%%
Size=25;
Downsampling_Ratio=100;


% %%
% NN=6;
% 
% plot(real(Data_Image_FFT(:,NN)));
% 

%


Mean_Array_Map_1D=Mean_Array_Map(:).*Mask_1D;
STD_Array_Map_1D=STD_Array_Map(:).*Mask_1D;
Mean_Array_Map_1D(isnan(Mean_Array_Map_1D))=[];
STD_Array_Map_1D(isnan(STD_Array_Map_1D))=[];
Mean_Array_Map_1D=downsample(Mean_Array_Map_1D,Downsampling_Ratio);
STD_Array_Map_1D=downsample(STD_Array_Map_1D,Downsampling_Ratio);

VAR_Array_Map_1D=STD_Array_Map_1D.^2;
scatter(Mean_Array_Map_1D,VAR_Array_Map_1D,Size,'filled');
xlabel('Signal (DN)','fontsize',15);
ylabel('Variance (DN^2)','fontsize',15);
set(gca,'fontsize',15)
%ylim([0 200]);
xlim([0 256]);
saveas(gcf,[Data_Save_Folder last_folder_name sprintf('_Method_%g',Noise_Removal_Method) '.png']);
