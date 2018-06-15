clear all
%%
root_folder_path='D:\170509_OctoPlus_PTC test\Diff Current\';
Data_Name='170509_OctoPlus_PTC test_DiffCurrent';
Current_Array=[150 200 250 300 400 500 600 700 800 900 1000];

Horizontal_Size=2048;
Size=25;

H_Offset=1;
H_Range=2047;
Horizontal_ROI=[H_Offset H_Offset+H_Range];


Mean_Array_Map=zeros(Horizontal_ROI(2)-Horizontal_ROI(1)+1,length(Current_Array));
STD_Array_Map=zeros(Horizontal_ROI(2)-Horizontal_ROI(1)+1,length(Current_Array));

for p=1:length(Current_Array)
    file_name=sprintf('%g',Current_Array(p));
    folder_path=[root_folder_path '\'];
    file_path=[folder_path file_name];

    Data_Save_Folder='D:\Processed Data\';
    

    Processed_Data_Path=[Data_Save_Folder Data_Name '_' file_name '.bin'];

    cd(folder_path);

    Data_Image=double(imread(file_path,'tiff'));
    
    Mean_Array_Map(Horizontal_ROI(1):Horizontal_ROI(2),p)=mean(Data_Image(:,Horizontal_ROI(1):Horizontal_ROI(2)),1);
    STD_Array_Map(Horizontal_ROI(1):Horizontal_ROI(2),p)=std(Data_Image(:,Horizontal_ROI(1):Horizontal_ROI(2)),1,1);
end

Mean_Array_Map_1D=Mean_Array_Map(:);
STD_Array_Map_1D=STD_Array_Map(:);
VAR_Array_Map_1D=STD_Array_Map_1D.^2;
scatter(Mean_Array_Map_1D,VAR_Array_Map_1D,Size,'filled');
xlabel('Signal (DN)','fontsize',15);
ylabel('Variance (DN^2)','fontsize',15);
set(gca,'fontsize',15)