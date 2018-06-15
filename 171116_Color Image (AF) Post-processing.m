clear all
%
%
folder_path='D:\171116_Post-processing of Scheme A Camera';
file_name='Image0531';
file_path=[folder_path '\' file_name '.tiff'];
new_file_path=[folder_path '\' file_name '_modi.tiff'];

Image_Original=imread(file_path,'tiff');
Image_Original=Image_Original(:,:,1:3);
imagesc(Image_Original);
axis equal
axis off

%
Binning_Factor=4;
Binned_Image=TSBinning(TSBinning(Image_Original,2,Binning_Factor),1,Binning_Factor);
imagesc(uint8(Binned_Image(:,:,1:3)));
axis equal
axis off
%

%
DC=80;
R_Min=0;
G_Min=0;
B_Min=0;

Ratio=1.9;
R_M=1;
G_M=1;
B_M=1;

Sub_Image=Binned_Image;
Sub_Image(:,:,1)=Ratio*R_M*(Sub_Image(:,:,1)-DC-R_Min);
Sub_Image(:,:,2)=Ratio*G_M*(Sub_Image(:,:,2)-DC-G_Min);
Sub_Image(:,:,3)=Ratio*G_M*(Sub_Image(:,:,3)-DC-B_Min);

subplot(2,1,1)
imagesc(Image_Original);
axis equal
axis off
subplot(2,1,2)
imagesc(uint8(Sub_Image));
axis equal
axis off

imwrite(uint8(Sub_Image),new_file_path);
