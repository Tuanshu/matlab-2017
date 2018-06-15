clear all
%%
Folder_Path='C:\TuanShu\171229_New FF DF\Test_1';
Row=1224;
Colomn=1024;

File_list=dir(Folder_Path);
Original_Length_file_list=length(File_list);
for p=1:Original_Length_file_list
    [pathstr,name,ext] = fileparts(File_list(Original_Length_file_list+1-p).name);
    if File_list(Original_Length_file_list+1-p).isdir ~= 0
        File_list(Original_Length_file_list+1-p)=[];
    elseif strcmp(File_list(Original_Length_file_list+1-p).name,'.') == 1
        File_list(Original_Length_file_list+1-p)=[];
    elseif strcmp(File_list(Original_Length_file_list+1-p).name,'..') == 1
        File_list(Original_Length_file_list+1-p)=[];
    elseif strcmp(ext,'.PNG')
        File_list(Original_Length_file_list+1-p)=[];
    elseif strcmp(File_list(Original_Length_file_list+1-p).name,'Processed Data') == 1
        File_list(Original_Length_file_list+1-p)=[];
    end
end

Frame=length(File_list);

   %% File Loading
Image_R=zeros(Colomn,Row);
Image_G=zeros(Colomn,Row);
Image_B=zeros(Colomn,Row);
for p=1:length(File_list)
    Image_Temp=imread([Folder_Path '\' File_list(p).name],'tif');
    Image_R=Image_R+double(Image_Temp(:,:,1));
    Image_G=Image_G+double(Image_Temp(:,:,2));
    Image_B=Image_B+double(Image_Temp(:,:,3));
    disp(p);
end
Image_R=Image_R/length(File_list);
Image_G=Image_G/length(File_list);
Image_B=Image_B/length(File_list);

Image_Flat_Field=zeros(Colomn,Row,3);

Image_Flat_Field(:,:,1)=Image_R;
Image_Flat_Field(:,:,2)=Image_G;
Image_Flat_Field(:,:,3)=Image_B;
Image_Flat_Field=uint8(Image_Flat_Field);
imagesc(Image_Flat_Field);
imwrite(Image_Flat_Field,[Folder_Path '\FF.tiff'],'tiff');
imwrite(Image_Flat_Field/2,[Folder_Path '\DF2.tiff'],'tiff');
imwrite(Image_Flat_Field/4,[Folder_Path '\DF4.tiff'],'tiff');
