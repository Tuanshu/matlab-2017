clear all

%%

Alpha_Image_Path='C:\TuanShu\170705_Raw Data Examine\Processed Data\20170704164625_0-290_gel_more_Npoint.tif';

Red_Image_Path='C:\TuanShu\170705_Raw Data Examine\Processed Data\20170704164625_0-290_gel_more_Npoint_Filtered_BW0.001_Offset-0.012.tif';
Green_Image_Path='C:\TuanShu\170705_Raw Data Examine\Processed Data\20170704164625_0-290_gel_more_Npoint_Filtered_BW0.003_Offset0.tif';
Blue_Image_Path='C:\TuanShu\170705_Raw Data Examine\Processed Data\20170704164625_0-290_gel_more_Npoint_Filtered_BW0.001_Offset0.012.tif';

Alpha=double(imread(Alpha_Image_Path));
RGB(:,:,1)=double(imread(Red_Image_Path));
RGB(:,:,3)=double(imread(Blue_Image_Path));
RGB(:,:,2)=RGB(:,:,3);%double(imread(Green_Image_Path));

Sum=sum(RGB,3);

RGB_Norm=RGB./repmat(Sum,[1 1 3]);

Alpha_RGB=repmat(Alpha/max(Alpha(:)),[1 1 3]).*RGB_Norm*2;

imagesc(Alpha_RGB);