clearvars 
%%
cd('C:\TuanShu\MATLAB\TSLib');
%%

Number=[1];%9
Set='Optronis';
%% PTC related param
Noise_Removal_Method=1; %1: do nothing, 2: normalization, 3: filtering, 4: 4-point

% Feature for Filtering mode
Filtering_Ratio=0.02;       %意思是, 若Array長度為A, 則將外側(低頻)之A*Filtering_Ratio個pixels設為0
% Feature for 4-point mode
Averaging_Factor=1;
N=4;          %for N-point
Unbias_Estimator=0.9213;

if strcmp(Set,'Optronis')
    root_folder_path='C:\TuanShu\';
    last_folder_name='171023_Euresys PTC Test (Optronis Camera)';
    If_Read_Header=1;
    Row=1280;
    Colomn=860;
    HeaderSize=8192;    %Byte
    Gap=5120;
    Bit_per_Pixel=8;
    Frame_Size=Row*Colomn*Bit_per_Pixel/8+Gap;
    Y_Offset=0;
end

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
    %% Header Reading
    if If_Read_Header == 1
        fileName=[folder_path '\' file_list_ori(1).name];   %Read the header if the 1st file
        headerInfo=NorpixHeaderReading(fileName)
        Row=headerInfo.ImageWidth;
        Colomn=headerInfo.ImageHeight;
        HeaderSize=headerInfo.HeaderSize;    %Byte
        Gap=headerInfo.TrueImageSize-headerInfo.ImageSizeBytes;
        Bit_per_Pixel=headerInfo.ImageBitDepth;
        Frame_Size=Row*Colomn*Bit_per_Pixel/8+Gap;
        Y_Offset=0;
        Frame_per_File_1stHeader=headerInfo.AllocatedFrames;    %just for reference
    end
    %% To Generate the real File List and Frame_Number_in_File_List Based on Data Size Estimation
    file_list=[];
    Frame_Number_in_File_List=[];
    for p=1:length(file_list_ori)
        Frame_Number_Calculated_Based_on_Size=(file_list_ori(p).bytes-HeaderSize)/Frame_Size;
        file_list=[file_list;repmat(file_list_ori(p),[Frame_Number_Calculated_Based_on_Size 1])]; 
        Frame_Number_in_File_List=[Frame_Number_in_File_List;[1:Frame_Number_Calculated_Based_on_Size]']; 
    end

    %%
    Frame=length(file_list);
    %%
    X=[1:Frame];
    Image_Stack=zeros(Row,Colomn,Frame);
    for p=1:Frame
        file_path=[folder_path file_list(p).name];
        fin=fopen(file_path);
        fseek(fin, HeaderSize+Frame_Size*(Frame_Number_in_File_List(p)-1), 'bof');
        %QQ=fread(fin,[Row,Colomn],'uint16');
        if  Bit_per_Pixel == 8
            Image_Stack(:,:,p)=fread(fin,[Row,Colomn],'uint8');
        elseif Bit_per_Pixel == 16
             Image_Stack(:,:,p)=fread(fin,[Row,Colomn],'uint16');
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
    All_Mean=mean(Mean_Array_for_norm);
    Mean_Stack_for_norm=repmat(Mean_Array_for_norm,[size(Image_Stack,1) size(Image_Stack,2) 1]);
    Image_Stack_norm=Image_Stack./Mean_Stack_for_norm*All_Mean;
    
    Image_Stack_FFT=fft((Image_Stack-Mean_Stack_for_norm),[],3);
    
    Image_Stack_FFT_Filter=Image_Stack_FFT;
    Image_Stack_FFT_Filter(:,:,1:Number_of_Filter_Pixel_one_side)=0;
    Image_Stack_FFT_Filter(:,:,(size(Image_Stack_FFT_Filter,1)-Number_of_Filter_Pixel_one_side+1):size(Image_Stack_FFT_Filter,1))=0;
    Image_Stack_Filter=real(ifft(Image_Stack_FFT_Filter,[],3))+Mean_Stack_for_norm;

    Averaged_Length=floor(size(Image_Stack,3)/Averaging_Factor);
    Temp=0;
    for q=1:Averaging_Factor
       Temp=Temp+Image_Stack(:,:,(Averaging_Factor-(q-1)):Averaging_Factor:(Averaging_Factor*Averaged_Length)-(q-1));
    end
    Image_Stack_Ave=Temp/Averaging_Factor;
    Image_Stack_Npoint=zeros(size(Image_Stack_Ave,1),size(Image_Stack_Ave,2),floor(size(Image_Stack_Ave,3)/N));
    for q=1:size(Image_Stack_Npoint,3)
        Temp=Image_Stack_Ave(:,:,((q-1)*N+1):(q*N));
        Image_Stack_Npoint(:,:,q)=((N*sum(Temp.^2,3)-sum(Temp,3).^2).^0.5)*(2^0.5)/N;
    end
    
    
    if Noise_Removal_Method == 1
        Mean_Array_Map(:,:,p)=mean(Image_Stack(X_ROI(1):X_ROI(2),Y_ROI(1):Y_ROI(2),:),3);
        STD_Array_Map(:,:,p)=std(Image_Stack(X_ROI(1):X_ROI(2),Y_ROI(1):Y_ROI(2),:),0,3);
    elseif Noise_Removal_Method == 2
        Mean_Array_Map(:,:,p)=mean(Image_Stack_norm(X_ROI(1):X_ROI(2),Y_ROI(1):Y_ROI(2),:),3);
        STD_Array_Map(:,:,p)=std(Image_Stack_norm(X_ROI(1):X_ROI(2),Y_ROI(1):Y_ROI(2),:),0,3);    
    elseif Noise_Removal_Method == 3
        Mean_Array_Map(:,:,p)=mean(Image_Stack_Filter(X_ROI(1):X_ROI(2),Y_ROI(1):Y_ROI(2),:),3);
        STD_Array_Map(:,:,p)=std(Image_Stack_Filter(X_ROI(1):X_ROI(2),Y_ROI(1):Y_ROI(2),:),0,3);        
    elseif Noise_Removal_Method == 4
        Mean_Array_Map(:,:,p)=mean(Image_Stack(X_ROI(1):X_ROI(2),Y_ROI(1):Y_ROI(2),:),3);
        STD_Array_Map(:,:,p)=mean(Image_Stack_Npoint(X_ROI(1):X_ROI(2),Y_ROI(1):Y_ROI(2),:),3)*Unbias_Estimator;        
    end
    disp(p);
    
end
 

%%
fclose('all');






clear all
%%
root_folder_path='C:\TuanShu\171023_Euresys PTC Test (Optronis Camera)\';
Current_Array=[0.1:0.1:0.7];


Horizontal_Size=2048;
Size=25;

H_Offset=300;
H_Range=30;
X_ROI=[H_Offset H_Offset+H_Range];


V_Offset=1;
V_Range=7;
Y_ROI=[V_Offset V_Offset+V_Range];

Noise_Removal_Method=4; %1: do nothing, 2: normalization, 3: filtering, 4: 4-point

% Feature for Filtering mode
Filtering_Ratio=0.02;       %意思是, 若Array長度為A, 則將外側(低頻)之A*Filtering_Ratio個pixels設為0
% Feature for 4-point mode
Averaging_Factor=1;
N=4;          %for N-point
Unbias_Estimator=0.9213;

Mean_Array_Map=zeros(X_ROI(2)-X_ROI(1)+1,Y_ROI(2)-Y_ROI(1)+1,length(Current_Array));
STD_Array_Map=zeros(X_ROI(2)-X_ROI(1)+1,Y_ROI(2)-Y_ROI(1)+1,length(Current_Array));

Mean_Array_for_norm=zeros(1,length(Current_Array));

Row=648;
Colomn=8;

for p=1:length(Current_Array)

   
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
ylim([0 200]);
xlim([0 4096]);
saveas(gcf,[Data_Save_Folder Data_Name sprintf('_Method_%g',Noise_Removal_Method) '.png']);
