clear all
%%
root_folder_path='D:\170328_Resolution test';

%root_folder_path='I:\161223_e2v EM1\CeYAG\Different Exp Time';

sub_folder_path='5.0 Pair_20x_no immersion_no Mirau_with FF';

sub_FF_path='FF';


Number=14;

folder_path=[root_folder_path '\' sub_folder_path '\'];

FF_path=[root_folder_path '\' sub_FF_path '\'];

Data_Name=sub_folder_path;
Display_Gain=4;


file_list=dir(folder_path);
file_list=file_list(3:end);
Stitched_Image=[];

FF_list=dir(FF_path);
FF_list=FF_list(3:end);
Stitched_Image=[];
%% FF
FF_name=FF_list(1).name;
FF_Raw=mean(double(imread([FF_path FF_name],'tiff')),3);

FF=1./FF_Raw;
FF=FF./max(FF(:));
%%

p=6;
N=1201;

file_name=file_list(p).name;
file_path=[folder_path file_name];

Data_Save_Folder='D:\170328_Resolution test\Processed Data\';
mkdir(Data_Save_Folder);
Processed_Data_Path=[Data_Save_Folder Data_Name '_' file_name '.bin'];
cd(folder_path);

Data_Image_Raw=mean(double(imread(file_path,'tiff')),3).*FF;
subplot(2,1,1);
imagesc(Data_Image_Raw)
axis equal
colormap('gray')
hold on
plot([1,size(Data_Image_Raw,2)],[N,N],'Color','r','LineWidth',2)

hold off
ylim([N-50 N+50])
xlim([400 550])

subplot(2,1,2);
plot(Data_Image_Raw(N,:),'-*');


%%
for p=1:length(file_list)
    file_name=file_list(p).name;
    file_path=[folder_path file_name];

    Data_Save_Folder='D:\170328_Resolution test\Processed Data\';
    mkdir(Data_Save_Folder);
    Processed_Data_Path=[Data_Save_Folder Data_Name '_' file_name '.bin'];
    cd(folder_path);

    Data_Image_Raw=mean(double(imread(file_path,'tiff')),3);
    imagesc(Data_Image_Raw)
end