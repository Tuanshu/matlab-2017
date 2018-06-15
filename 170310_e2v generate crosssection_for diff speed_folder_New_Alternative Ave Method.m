clear all
%%
root_folder_path='F:\170310_EM1.4 test\Chnage Anamorphic\Religned';
sub_folder_path='Forearm';
Number=15;

%root_folder_path='I:\161223_e2v EM1\CeYAG\Different Exp Time';

Averaging_Factor=16;
Vertical_Maximum_Depth=1400;    %pixel 1285

folder_path=[root_folder_path '\' sub_folder_path '_' sprintf('%g',Number) '\'];


Data_Name=sub_folder_path;
it=29;
Gain=-24;
Display_Gain=4;
Horizontal_Size=1024;
Size=25;


SR_Lateral=0.84;
SR_Axial=0.279;
Axial_ave_Factor=1;
%

H_Offset_View=130;
H_Range_View=650;
Horizontal_ROI_View=[H_Offset_View H_Offset_View+H_Range_View];

%
H_Offset=360;
H_Range=360;
Horizontal_ROI=[H_Offset H_Offset+H_Range];

If_norm=1; %0: do nothing, 1: normalization

% Feature for Filtering mode
Filtering_Ratio=0.02;       %意思是, 若Array長度為A, 則將外側(低頻)之A*Filtering_Ratio個pixels設為0
% Feature for 4-point mode
Product_of_Axial_Decimation_Factor_and_Ave_Factor=Averaging_Factor;
N=4;          %for N-point
Unbias_Estimator=0.9213;
% 
% Mean_Array_Map=zeros(Horizontal_ROI(2)-Horizontal_ROI(1)+1,length(Current_Array));
% STD_Array_Map=zeros(Horizontal_ROI(2)-Horizontal_ROI(1)+1,length(Current_Array));
% 
% Mean_Array_for_norm=zeros(1,length(Current_Array));



Current_Array=[1];%[3:3:24];

file_list=dir(folder_path);
file_list=file_list(4:end);
Stitched_Image=[];
for p=1:length(file_list)
    file_name=file_list(p).name;
    file_path=[folder_path '\' file_name];

    Data_Save_Folder='F:\Processed Data\';

    Processed_Data_Path=[Data_Save_Folder Data_Name '_' file_name '.bin'];

    cd(folder_path);

    Data_Image_Raw=double(imread(file_path,'tiff'));
    Data_Image_Raw_Decimated=Data_Image_Raw(1:floor(Product_of_Axial_Decimation_Factor_and_Ave_Factor/Averaging_Factor):end);
       
    if If_norm == 0
           Data_Image=Data_Image_Raw;
    elseif If_norm == 1
        Mean_Array_for_norm=mean(Data_Image_Raw(:,Horizontal_ROI(1):Horizontal_ROI(2)),2); % 改成可調ROI
        All_Mean=mean(Mean_Array_for_norm);
        Mean_Image_for_norm=repmat(Mean_Array_for_norm,[1 size(Data_Image_Raw,2)]);
        Data_Image=Data_Image_Raw./Mean_Image_for_norm*All_Mean;
    end
    
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
    Stitched_Image=[Stitched_Image; Data_Image_Npoint];
end

%imagesc(Data_Image_Npoint);
%
if size(Stitched_Image,1)>Vertical_Maximum_Depth
    Stitched_Image=Stitched_Image(1:Vertical_Maximum_Depth,:);
end

Temp=0;
Axial_Length_Original=size(Stitched_Image,1);
Axial_Length_Used=floor(Axial_Length_Original/Axial_ave_Factor)*Axial_ave_Factor;
Reduced_Length=Axial_Length_Used/Axial_ave_Factor;
for p=1:Axial_ave_Factor
   Temp=Temp+Stitched_Image((Axial_ave_Factor-(p-1)):Axial_ave_Factor:(Axial_ave_Factor*Reduced_Length)-(p-1),:);
end
Data_Image_Reduced=Temp/Axial_ave_Factor;




%
C_max=24;   %32
C_min=1.6;  %0.8

Data_Image_Npoint_normalized=Display_Gain*(Data_Image_Reduced(:,Horizontal_ROI_View(1):Horizontal_ROI_View(2))-C_min)/(C_max-C_min);
Data_Image_Npoint_normalized(Data_Image_Npoint_normalized<0)=0;
Data_Image_Npoint_normalized(Data_Image_Npoint_normalized>1)=1;



subplot(1,1,1)
imagesc(Data_Image_Npoint_normalized);
%axis off
%axis equal
caxis([0 1]);
X_Label_Array=[0  100 200 300 400 500];
X_Label_Index_Array=X_Label_Array/SR_Lateral;

Y_Label_Array=[0 50 100 150 200 250 300 350 400];
Y_Label_Index_Array=round(Y_Label_Array/SR_Axial/Axial_ave_Factor);
set(gca,'xtick',X_Label_Index_Array,'xticklabel',X_Label_Array,'ytick',Y_Label_Index_Array,'yticklabel',Y_Label_Array);
xlabel('Lateral Position (micron)');
ylabel('Axial Position (micron)');
colormap(gray);




if If_norm ==0
    imwrite(Data_Image_Npoint_normalized,[Processed_Data_Path sprintf('_Tissue %g',Axial_ave_Factor) sprintf('_Ave %g',Number) '_Bscan.png'],'png');
elseif If_norm ==1
    imwrite(Data_Image_Npoint_normalized,[Processed_Data_Path sprintf('_Tissue %g',Axial_ave_Factor) sprintf('_Ave %g',Number) '_Bscan_norm.png'],'png');
end
    


% %%
% NN=6;
% 
% plot(real(Data_Image_FFT(:,NN)));
% 

%

% 
% Mean_Array_Map_1D=Mean_Array_Map(:);
% STD_Array_Map_1D=STD_Array_Map(:);
% VAR_Array_Map_1D=STD_Array_Map_1D.^2;
% scatter(Mean_Array_Map_1D,VAR_Array_Map_1D,Size,'filled');
% xlabel('Signal (DN)','fontsize',15);
% ylabel('Variance (DN^2)','fontsize',15);
% set(gca,'fontsize',15)
% ylim([0 100]);
% xlim([0 4000]);