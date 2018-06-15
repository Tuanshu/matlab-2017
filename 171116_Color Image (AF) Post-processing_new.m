clear all
%
%
folder_path='D:\171116_Post-processing of Scheme A Camera\long exposure\0.5it';
file_name='Image0738';
file_path=[folder_path '\' file_name '.tiff'];
new_file_path=[folder_path '\' file_name '_modi.tiff'];

Image_Original=imread(file_path,'tiff');
Image_Original=Image_Original(:,:,1:3);
% imagesc(Image_Original);
% axis equal
% axis off

%% Separate Intensity and Color
Image_Intensity=mean(Image_Original,3);
Image_NormalizedColor=double(Image_Original)./repmat(Image_Intensity,[1 1 3])/3;

%% Filtering
Window_Size=9;
Sigma=2;


Window_Size_Color=9;
Sigmae_Color=2;

GuassFilter = fspecial('gaussian',Window_Size,Sigma);
GuassFilter_Color = fspecial('gaussian',Window_Size_Color,Sigmae_Color);

Image_Intensity_Filtered=filter2(GuassFilter,Image_Intensity,'same');

for p=1:3
    Image_NormalizedColor_Filtered(:,:,p)=filter2(GuassFilter_Color,Image_NormalizedColor(:,:,p),'same');

end
%
Image_Reconstructed_1=uint8(Image_NormalizedColor_Filtered.*repmat(Image_Intensity_Filtered,[1 1 3])*3);

% imagesc(Image_Reconstructed_1);
% axis equal
% axis off
%% Mono Contrast Enhacement
DC=30;
Ratio=1.2;
Image_Intensity_Modified=Ratio.*(Image_Intensity_Filtered-DC);
Image_Reconstructed_2=uint8(Image_NormalizedColor_Filtered.*repmat(Image_Intensity_Modified,[1 1 3])*3);

% subplot(2,1,1)
% imagesc(Image_Original);
% xlim([1 size(Image_Original,1)])
% ylim([1 size(Image_Original,2)])
% axis equal
% axis off
% subplot(2,1,2)
% imagesc(uint8(Image_Reconstructed_2));
% xlim([1 size(Image_Reconstructed_2,1)])
% ylim([1 size(Image_Reconstructed_2,2)])
% axis equal
% axis off

%% Color Tunning
% Reference_Path='D:\171116_Dermoscope Reference\IMG_0051_Test.jpg';
% RGB_modi=ColorMimicking(Image_Reconstructed_2,Reference_Path,'tiff','jpg');
R_Ratio=1.1975;%RGB_modi(1);
G_Ratio=0.9665;%RGB_modi(2);
B_Ratio=0.8354;%RGB_modi(3);
ColorSpace=zeros(1,1,3);
ColorSpace(1)=R_Ratio/(R_Ratio+G_Ratio+B_Ratio)*3;
ColorSpace(2)=G_Ratio/(R_Ratio+G_Ratio+B_Ratio)*3;
ColorSpace(3)=B_Ratio/(R_Ratio+G_Ratio+B_Ratio)*3;


Image_NormalizedColor_Tunned=Image_NormalizedColor_Filtered.*repmat(ColorSpace,[size(Image_NormalizedColor_Filtered,1) size(Image_NormalizedColor_Filtered,2) 1]);

Image_Reconstructed_3=uint8(Image_NormalizedColor_Tunned.*repmat(Image_Intensity_Modified,[1 1 3])*3);
% 
% subplot(2,1,1)
% imagesc(Image_Original);
% xlim([1 size(Image_Original,1)])
% ylim([1 size(Image_Original,2)])
% axis equal
% axis off
% subplot(2,1,2)
% imagesc(uint8(Image_Reconstructed_3));
% xlim([1 size(Image_Reconstructed_3,1)])
% ylim([1 size(Image_Reconstructed_3,2)])
% axis equal
% axis off

%%
%imwrite(uint8(Image_Reconstructed_3),new_file_path);
