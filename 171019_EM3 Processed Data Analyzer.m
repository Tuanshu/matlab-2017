clear all
%%

Data_Folder_Path='C:\TuanShu\171018_EM2 E-scan V Test\Processed Data\';
File_Name='20171018180152_1x_3.5it_V_Back to 8.3mW_3D.raw';
Row=1024;
Colomn=1024;

X_Range=[313 712];%[155 868];
Y_Range=[313 712];%[155 868];
Z_Range=[101 500];%[1 714];
%%
file_path=[Data_Folder_Path File_Name];
File_Info=dir(file_path);

Frame=File_Info.bytes*8/64/Row/Colomn;  %64 for double
%%
Raw_Data=zeros(Row,Colomn,Frame);
fin=fopen(file_path);
for p=1:Frame
    %fseek(fin, Byte_Skip+Row*Colomn*2*(Frame_Number_in_File_List(Current_Count)-1), 'bof');
    Raw_Data(:,:,p)=fread(fin,[Row,Colomn],'double'); %*Frame   不知為何, 看起來就像是要除16
    disp(p);
end
fclose(fin);
%imagesc(Raw_Data(:,:,1));
%% Sub Volume
Raw_Data_Cut=Raw_Data(X_Range(1):X_Range(2),Y_Range(1):Y_Range(2),Z_Range(1):Z_Range(2));
clear Raw_Data
%%
imagesc(squeeze(Raw_Data_Cut(:,370,:)))
%% Contrast Analysis (refer to Contrast Measure based on Squared Laplacian (CMSL) 
Gx_Array=zeros([Frame 1]);
Gy_Array=zeros([Frame 1]);
for p=1:(Z_Range(2)-Z_Range(1)+1)
    Gx_Map=abs(diff(Raw_Data_Cut(:,:,p),1,1));
    Gy_Map=abs(diff(Raw_Data_Cut(:,:,p),1,2));
    Gx_Array(p)=sum(Gx_Map(:).^2);
    Gy_Array(p)=sum(Gy_Map(:).^2);
    disp(p);
end

plot(1:Frame,Gx_Array,1:Frame,Gy_Array)

%% Spatial Frequency Analysis (XY)
FFT_Map_XY=zeros(X_Range(2)-X_Range(1)+1,Y_Range(2)-Y_Range(1)+1);
for p=1:(Z_Range(2)-Z_Range(1)+1)
    Temp=abs(fftshift(fftshift(fft2(Raw_Data_Cut(:,:,p)),1),2));
    Temp=(Temp-min(Temp(:)))/(max(Temp(:))-min(Temp(:)));
    FFT_Map_XY=FFT_Map_XY+Temp;
    disp(p);
end
FFT_Map_XY=FFT_Map_XY/Frame;
%%
imagesc(FFT_Map_XY)
axis equal
caxis([0 0.001])
%% Spatial Frequency Analysis (XZ)
FFT_Map_XZ=zeros(X_Range(2)-X_Range(1)+1,Z_Range(2)-Z_Range(1)+1);
for p=1:(Y_Range(2)-Y_Range(1)+1)
    Temp=abs(fftshift(fftshift(fft2(squeeze(Raw_Data_Cut(:,p,:))),1),2));
    Temp=(Temp-min(Temp(:)))/(max(Temp(:))-min(Temp(:)));
    FFT_Map_XZ=FFT_Map_XZ+Temp;
    disp(p);
end
FFT_Map_XZ=FFT_Map_XZ/Frame;
%%
imagesc(FFT_Map_XZ)
axis equal
caxis([0 0.01])