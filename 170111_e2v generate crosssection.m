clear all
%%
root_folder_path='F:\170111_e2v EM1 new setup\Different Gain';
%root_folder_path='I:\161223_e2v EM1\CeYAG\Different Exp Time';

Data_Name='170111_e2v EM1 new setup';
Current_Array=[84];%[3:3:24];
it=29;
Gain=-24;
Display_Gain=4;
Horizontal_Size=1024;
Size=25;


SR_Lateral=0.888;
SR_Axial=0.279;
Axial_ave_Factor=3;

%

H_Offset_View=200;
H_Range_View=700;
Horizontal_ROI_View=[H_Offset_View H_Offset_View+H_Range_View];

%
H_Offset=360;
H_Range=460;
Horizontal_ROI=[H_Offset H_Offset+H_Range];

If_norm=1; %1: do nothing, 2: normalization, 3: filtering, 4: 4-point

% Feature for Filtering mode
Filtering_Ratio=0.02;       %�N��O, �YArray���׬�A, �h�N�~��(�C�W)��A*Filtering_Ratio��pixels�]��0
% Feature for 4-point mode
Averaging_Factor=16;
Product_of_Axial_Decimation_Factor_and_Ave_Factor=16;
N=4;          %for N-point
Unbias_Estimator=0.9213;

Mean_Array_Map=zeros(Horizontal_ROI(2)-Horizontal_ROI(1)+1,length(Current_Array));
STD_Array_Map=zeros(Horizontal_ROI(2)-Horizontal_ROI(1)+1,length(Current_Array));

Mean_Array_for_norm=zeros(1,length(Current_Array));

for p=1:length(Current_Array)
    file_name=sprintf('Tissue_%g_%gdB_it%g',Current_Array(p),Gain,it);
    folder_path=[root_folder_path '\'];
    file_path=[folder_path file_name];

    Data_Save_Folder='F:\Processed Data\';

    Processed_Data_Path=[Data_Save_Folder Data_Name '_' file_name '.bin'];

    cd(folder_path);

    Data_Image_Raw=double(imread(file_path,'tiff'));
    Data_Image_Raw_Decimated=Data_Image_Raw(1:floor(Product_of_Axial_Decimation_Factor_and_Ave_Factor/Averaging_Factor):end);
       
    if If_norm == 0
        Mean_Array_Map(Horizontal_ROI(1):Horizontal_ROI(2),p)=mean(Data_Image_Raw(:,Horizontal_ROI(1):Horizontal_ROI(2)),1);
        STD_Array_Map(Horizontal_ROI(1):Horizontal_ROI(2),p)=std(Data_Image_Raw(:,Horizontal_ROI(1):Horizontal_ROI(2)),1,1);
        Data_Image=Data_Image_Raw;
    elseif If_norm == 1
        Mean_Array_for_norm=mean(Data_Image_Raw(:,Horizontal_ROI(1):Horizontal_ROI(2)),2); % �令�i��ROI
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
    
end

%imagesc(Data_Image_Npoint);
%
Temp=0;
Axial_Length_Original=size(Data_Image_Npoint,1);
Axial_Length_Used=floor(Axial_Length_Original/Axial_ave_Factor)*Axial_ave_Factor;
Reduced_Length=Axial_Length_Used/Axial_ave_Factor;
for p=1:Axial_ave_Factor
   Temp=Temp+Data_Image_Npoint((Axial_ave_Factor-(p-1)):Axial_ave_Factor:(Axial_ave_Factor*Reduced_Length)-(p-1),:);
end
Data_Image_Reduced=Temp/Axial_ave_Factor;



%
C_max=32;
C_min=0.8;

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

Y_Label_Array=[0 50 100 150 200 250 300 350];
Y_Label_Index_Array=round(Y_Label_Array/SR_Axial/Axial_ave_Factor);
set(gca,'xtick',X_Label_Index_Array,'xticklabel',X_Label_Array,'ytick',Y_Label_Index_Array,'yticklabel',Y_Label_Array);
xlabel('Lateral Position (micron)');
ylabel('Axial Position (micron)');
colormap(gray);
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