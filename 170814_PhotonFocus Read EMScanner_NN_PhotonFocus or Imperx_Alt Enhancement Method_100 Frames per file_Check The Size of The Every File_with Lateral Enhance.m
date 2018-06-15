clearvars 
%%

root_folder_path='C:\TuanShu\';
last_folder_name='170811_80PROJ_it0.32_030-050';
Number=[];

System='780';
Camera='MV1';    %MV1 or B0620;

If_EM2=1;
If_Auto_Brightness=1;
Auto_Brightness_TH=0.99;
If_FolderNaming=1;
N=4;

If_Calculate_NN=1;
If_Search_Glass=1;
If_Save_Raw=0;
Assumed_Glass_Interface_Index=711;
If_Norm=1;
If_Enhace_Deep_Signal=1;
If_Enhace_Lateral_Signal=1;
BI=-0.1;
BF=1;
% Data format related
Row=1024;
Colomn=11;
Y_Offset=0;
X_ROI=[1 1024];%[1 1024];%[300 700];%[1 1024];%[80 80+856-1];%[35 35+781-1];%[80 80+856-1];
NN_ROI=[300 700];
NN_Y_Range=4;
% 
if strcmp(Camera,'MV1')
    Number_of_Frame_per_File=1000;  %2
elseif strcmp(Camera,'B0620')
    Number_of_Frame_per_File=1;
end
Byte_Skip=Y_Offset*Row*2;
% Processing related
Column_Binning_Factor=1;
Row_Binning_Factor=1;

if strcmp(Camera,'MV1')
    Lateral_Ave_Factor=1;   %New Param 170413
elseif strcmp(Camera,'B0620')
    Lateral_Ave_Factor=2;   %New Param 170413
end


if strcmp(System,'560')
    Axial_ave_Factor=4;
elseif strcmp(System,'780')    
    if strcmp(Camera,'MV1')
        Axial_ave_Factor=2;%2;
    elseif strcmp(Camera,'B0620')
        Axial_ave_Factor=3;
    end
end

if strcmp(Camera,'MV1')
    Product_of_Axial_Decimation_Factor_and_Ave_Factor=2;%round(2);
elseif strcmp(Camera,'B0620')
    Product_of_Axial_Decimation_Factor_and_Ave_Factor=round(4);
end
ave_factor=Product_of_Axial_Decimation_Factor_and_Ave_Factor;%[4*4*4];

Phase_Shift=0;

NN_Window_Size=100;
Signal_Window_Size=20;
Separate_Size=25;

% For Last File Number Estimation
Bit_per_Pixel=16;
Frame_Size=Row*Colomn*Bit_per_Pixel/8;  %Byte


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
for QQQ=1:length(folder_list)
    folder_path=[parent_folder_path '\' folder_list(QQQ).name '\'];
    cd(folder_path);
    micron_per_frame=0.2/ave_factor/N;
    Offset_1=0;
    Offset_2=0; 

    Max_Number_of_Frame=3000000;

    %%
    file_list_ori=dir(folder_path);
    Original_Length_file_list=length(file_list_ori);

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
    
    %% To Generate the real File List and Frame_Number_in_File_List Based on Data Size Estimation
    file_list=[];
    Frame_Number_in_File_List=[];
    for p=1:length(file_list_ori)
        Frame_Number_Calculated_Based_on_Size=file_list_ori(p).bytes/Frame_Size;
        file_list=[file_list;repmat(file_list_ori(p),[Frame_Number_Calculated_Based_on_Size 1])]; 
        Frame_Number_in_File_List=[Frame_Number_in_File_List;[1:Frame_Number_Calculated_Based_on_Size]']; 
    end

    %%
    file_list=downsample(file_list,floor(Product_of_Axial_Decimation_Factor_and_Ave_Factor/ave_factor),Phase_Shift);
    Frame_Number_in_File_List=downsample(Frame_Number_in_File_List,floor(Product_of_Axial_Decimation_Factor_and_Ave_Factor/ave_factor),Phase_Shift);
    Frame=length(file_list);
    %%
    Ave_Temp=zeros(Row,Colomn,ave_factor);
    After_Npoint_Frame_Length=floor(Frame/N/ave_factor);
    After_Npoint_Image_Stack=zeros(Row/Row_Binning_Factor,Colomn/Column_Binning_Factor,min(Max_Number_of_Frame,After_Npoint_Frame_Length));
    if If_Calculate_NN ==1
        Intensity_Stack=zeros(Row/Row_Binning_Factor,Colomn/Column_Binning_Factor,min(Max_Number_of_Frame,After_Npoint_Frame_Length));
    end
    

    X=[1:Frame];
    Frame_Ave=zeros([min(Max_Number_of_Frame,After_Npoint_Frame_Length)*N 1]);
    for p=1:min(Max_Number_of_Frame,After_Npoint_Frame_Length)
        Npoint_Temp=zeros(Row/Row_Binning_Factor,Colomn/Column_Binning_Factor,N);
        for q=1:N
            for r=1:ave_factor
                %%      
                Current_Count=(p-1)*N*ave_factor+(q-1)*ave_factor+r;
                file_path=[folder_path file_list(Current_Count).name];

                fin=fopen(file_path);

                fseek(fin, Byte_Skip+Row*Colomn*2*(Frame_Number_in_File_List(Current_Count)-1), 'bof');
                %QQ=fread(fin,[Row,Colomn],'uint16');   
                Ave_Temp(:,:,r)=fread(fin,[Row,Colomn],'uint16'); %*Frame   不知為何, 看起來就像是要除16
                %%
                %fclose(fin);
                fclose('all');
            end
            Mean=mean(Ave_Temp,3);

            % Row Binning
            Temp=0;
            Reduced_Width=size(Mean,1)/Row_Binning_Factor;
            for tt=1:Row_Binning_Factor
                Temp=Temp+Mean((Row_Binning_Factor-(tt-1)):Row_Binning_Factor:(Row_Binning_Factor*Reduced_Width)-(tt-1),:);
            end
            Temp=Temp/Row_Binning_Factor;

            Reduced_Length=size(Mean,2)/Column_Binning_Factor;
            for tt=1:Column_Binning_Factor
                Npoint_Temp(:,:,q)=Npoint_Temp(:,:,q)+Temp(:,(Column_Binning_Factor-(tt-1)):Column_Binning_Factor:(Column_Binning_Factor*Reduced_Length)-(tt-1));
            end
            Npoint_Temp(:,:,q)=Npoint_Temp(:,:,q)/Column_Binning_Factor;


            Frame_Ave((p-1)*N+q)=mean(mean(Npoint_Temp(:,:,q)));    %Always calculated

            if If_Norm ==1
                Npoint_Temp(:,:,q)=Npoint_Temp(:,:,q)/Frame_Ave((p-1)*N+q);
            end
        end

        After_Npoint_Image_Stack(:,:,p)=((N*sum(Npoint_Temp.^2,3)-sum(Npoint_Temp,3).^2).^0.5)*(2^0.5)/N;

        if If_Calculate_NN ==1
            Intensity_Stack(:,:,p)=mean(Npoint_Temp,3);
        end

        disp(p);
    end

    if If_Norm ==1
        After_Npoint_Image_Stack=After_Npoint_Image_Stack*mean(Frame_Ave);
        if If_Calculate_NN ==1
            Intensity_Stack=Intensity_Stack*mean(Frame_Ave);
        end
    end
    
    if If_Calculate_NN ==1
        %% Glass interface detection
        if If_Search_Glass ==1
            [temp_value Glass_Interface_Index]=max(squeeze(mean(mean(After_Npoint_Image_Stack,1),2)));
            if Glass_Interface_Index+(NN_Window_Size+Separate_Size)>size(After_Npoint_Image_Stack,3)
                Glass_Interface_Index=Assumed_Glass_Interface_Index;
            end
        else
            Glass_Interface_Index=Assumed_Glass_Interface_Index;
        end
        if If_EM2 ==1
            NN_Range=[Glass_Interface_Index+Separate_Size+1 Glass_Interface_Index+Separate_Size+NN_Window_Size];
            Signal_Range=[Glass_Interface_Index-Signal_Window_Size/2+1 Glass_Interface_Index+Signal_Window_Size/2];
        else
            NN_Range=[Glass_Interface_Index-Separate_Size-NN_Window_Size Glass_Interface_Index-Separate_Size-1];
            Signal_Range=[Glass_Interface_Index-Signal_Window_Size/2+1 Glass_Interface_Index+Signal_Window_Size/2];
        end
        NN_Map=std(After_Npoint_Image_Stack(NN_ROI(1):NN_ROI(2),NN_Y_Range,NN_Range(1):NN_Range(2))./Intensity_Stack(NN_ROI(1):NN_ROI(2),NN_Y_Range,NN_Range(1):NN_Range(2)),0,3);
        IE_Map=max(After_Npoint_Image_Stack(NN_ROI(1):NN_ROI(2),NN_Y_Range,Signal_Range(1):Signal_Range(2))./Intensity_Stack(NN_ROI(1):NN_ROI(2),NN_Y_Range,Signal_Range(1):Signal_Range(2)),[],3);
        Signal_Map=max(After_Npoint_Image_Stack(NN_ROI(1):NN_ROI(2),NN_Y_Range,Signal_Range(1):Signal_Range(2)),[],3);
        Noise_Map=std(After_Npoint_Image_Stack(NN_ROI(1):NN_ROI(2),NN_Y_Range,NN_Range(1):NN_Range(2)),0,3);
        DF_Map=mean(After_Npoint_Image_Stack(NN_ROI(1):NN_ROI(2),NN_Y_Range,NN_Range(1):NN_Range(2)),3);
        Intensity_Map_Noise=mean(Intensity_Stack(NN_ROI(1):NN_ROI(2),NN_Y_Range,NN_Range(1):NN_Range(2)),3);
        Intensity_Map_Sig=mean(Intensity_Stack(NN_ROI(1):NN_ROI(2),NN_Y_Range,Signal_Range(1):Signal_Range(2)),3);
        NN=mean(NN_Map(:));
        IE=max(IE_Map(:));
        DF=mean(DF_Map(:));
        Signal=max(Signal_Map(:));
        Noise=mean(Noise_Map(:));
        Intensity=mean(Intensity_Map_Noise(:));
        IE_Sub_Map=((Signal_Map.^2-DF.^2).^0.5)./Intensity_Map_Sig;
        IE_Sub=max(IE_Sub_Map(:));
    end
    
    
    disp(QQQ);
    
    %%
    Maximum_Axial_Frame=20000;
    Temp=0;
    Axial_Length_Original=size(After_Npoint_Image_Stack,3);
    Axial_Length_Used=floor(Axial_Length_Original/Axial_ave_Factor)*Axial_ave_Factor;
    Reduced_Length=Axial_Length_Used/Axial_ave_Factor;
    for p=1:Axial_ave_Factor
       Temp=Temp+After_Npoint_Image_Stack(:,:,(Axial_ave_Factor-(p-1)):Axial_ave_Factor:(Axial_ave_Factor*Reduced_Length)-(p-1));
    end
    Reduced_Stack=Temp/Axial_ave_Factor;
    
    %Reduced_Image=squeeze(Reduced_Stack(:,size(Reduced_Stack,2)/2,:))';
    %Reduced_Image=squeeze(mean(Reduced_Stack(:,1:Y_Ave_Factor/Column_Binning_Factor,:),2))';
    Reduced_Image=squeeze(mean(Reduced_Stack(:,:,:),2))';
    disp(size(Reduced_Stack,2))
    Reduced_Image=Reduced_Image(1:min(size(Reduced_Image,1),Maximum_Axial_Frame),X_ROI(1):X_ROI(2));

    if If_EM2 ==1
        Reduced_Image=Reduced_Image(size(Reduced_Image,1):-1:1,:);
    end   
   
      %%
    if If_Enhace_Lateral_Signal
        XX=[1:size(Reduced_Image,2)]';
        FitResult=fit(XX,mean(Reduced_Image,1)','poly2');
        plot(XX,mean(Reduced_Image,1),XX,FitResult.p1*XX.^2 + FitResult.p2*XX + FitResult.p3);
        Lateral_Profile=FitResult.p1*XX.^2 + FitResult.p2*XX + FitResult.p3;
        
        
        Lateral_Enhace_Coef=mean(Lateral_Profile)./(Lateral_Profile);

        
        Reduced_Image=Reduced_Image.*repmat(Lateral_Enhace_Coef',[size(Reduced_Image,1) 1]);        
    end    
    
    
    %%
    
    if If_Enhace_Deep_Signal

        Linear_Enhace_Coef=[BI:(BF-BI)/(size(Reduced_Image,1)-1):BF]';

        
        Reduced_Image=Reduced_Image.*repmat(Linear_Enhace_Coef,[1 size(Reduced_Image,2)]);
        
    end    

    
%%        
    %clear After_Npoint_Image_Stack
    

     
    
    if If_Auto_Brightness == 1
        nbin=250;
        Histogram=hist(Reduced_Image(:),nbin);

        Total=sum(Histogram);
        Int_Hist = cumsum(Histogram)/Total;

        plot(Int_Hist);

        C_max_Auto=find(Int_Hist>Auto_Brightness_TH,1,'first')/nbin*max(Reduced_Image(:));
        C_max=C_max_Auto;
        C_min=C_max*0.05;
    else
    %% Hist for auto brightness



        C_max=30;
        C_min=1.5;
    
    end
    Reduced_Image_normalized=(Reduced_Image-C_min)/(C_max-C_min);
    Reduced_Image_normalized(Reduced_Image_normalized<0)=0;
    Reduced_Image_normalized(Reduced_Image_normalized>1)=1;
    
    
    %% Lateral Ave
    Temp=0;
    Lateral_Length_Original=size(Reduced_Image_normalized,2);
    Reduced_Length=Lateral_Length_Original/Lateral_Ave_Factor;
    for p=1:Lateral_Ave_Factor
       Temp=Temp+Reduced_Image_normalized(:,(Lateral_Ave_Factor-(p-1)):Lateral_Ave_Factor:(Lateral_Ave_Factor*Reduced_Length)-(p-1));
    end
    Lateral_Binned_Image=Temp/Lateral_Ave_Factor;
    
    %Reduced_Image=squeeze(Reduced_Stack(:,size(Reduced_Stack,2)/2,:))';
    %Reduced_Image=squeeze(mean(Reduced_Stack(:,1:Y_Ave_Factor,:),2))';
    
    
    
    
    
    %
    
    subplot(1,1,1)
    imagesc(Lateral_Binned_Image);
    axis off
    axis equal
    caxis([0 1]);
    set(gca,'xtick',[0  20  40  60  80  100],'xticklabel',[0  20  40  60  80  100]*0.33,'ytick',[0  100  200  300  400  500 600],'yticklabel',[0  100  200  300  400  500 600]*0.2645);
    xlabel('Lateral Position (micron)');
    ylabel('Axial Position (micron)');
    colormap(gray);
    if If_Calculate_NN ==1
        if Glass_Interface_Index+(NN_Window_Size+Separate_Size)>size(After_Npoint_Image_Stack,3)
            Test_Image=squeeze(After_Npoint_Image_Stack(X_ROI(1):X_ROI(2),NN_Y_Range,:))';

        else
            Test_Image=squeeze(After_Npoint_Image_Stack(X_ROI(1):X_ROI(2),NN_Y_Range,NN_Range(1):NN_Range(2)))';
        end
        Test_Image=Test_Image/20;
    end
    if If_Enhace_Deep_Signal
        Processed_Sub_Folder_Name=[Data_Save_Folder sprintf('System_%s_AAF_%d_LAF_%d_RBF_%d_CBF_%d_ABTH_%g_AVE_%d_Phase_%d_ProductAVE_%d_Colomn_%d_DeepEnhace_XWidth%d_BI%g_BF%g',System,Axial_ave_Factor,Lateral_Ave_Factor,Row_Binning_Factor,Column_Binning_Factor,Auto_Brightness_TH,ave_factor,Phase_Shift,Product_of_Axial_Decimation_Factor_and_Ave_Factor,Colomn,X_ROI(2)-X_ROI(1),BI,BF) '\'];
        if exist(Processed_Sub_Folder_Name)==0
            mkdir(Processed_Sub_Folder_Name);
        end
        if length(Number)>0
            if If_FolderNaming ==0
                imwrite(Lateral_Binned_Image,[Processed_Sub_Folder_Name sprintf('%03d_%s.png',Number(QQQ),System)],'png');
            else
                imwrite(Lateral_Binned_Image,[Processed_Sub_Folder_Name folder_list(QQQ).name '.png'],'png');
            end
            if If_Save_Raw == 1
                fid = fopen([Processed_Sub_Folder_Name sprintf('%03d_%s.raw',Number(QQQ),System)], 'w+');
                fwrite(fid, Lateral_Binned_Image', 'single');
                fclose(fid);
            end
            if If_Calculate_NN ==1
                dlmwrite([Processed_Sub_Folder_Name 'Inten_NN Record.txt'],[Number(QQQ) Intensity NN],'delimiter','\t','newline','pc','-append');
                dlmwrite([Processed_Sub_Folder_Name 'IE_Sub Record.txt'],[Number(QQQ) IE_Sub],'delimiter','\t','newline','pc','-append');
                imwrite(Test_Image,[Processed_Sub_Folder_Name sprintf('NN_%03d_%s.png',Number(QQQ),System)],'png');
            end
        elseif Y_Offset~=0
            if If_FolderNaming ==0
                imwrite(Lateral_Binned_Image,[Processed_Sub_Folder_Name sprintf('%03d_%s.png',QQQ,System)],'png');
            else
                imwrite(Lateral_Binned_Image,[Processed_Sub_Folder_Name folder_list(QQQ).name '.png'],'png');
            end
            if If_Save_Raw == 1
                fid = fopen([Processed_Sub_Folder_Name sprintf('%03d_%s.raw',QQQ,System)], 'w+');
                fwrite(fid, Lateral_Binned_Image', 'single');
                fclose(fid);
            end
            if If_Calculate_NN ==1
                dlmwrite([Processed_Sub_Folder_Name 'Inten_NN Record.txt'],[QQQ Intensity NN],'delimiter','\t','newline','pc','-append');
                dlmwrite([Processed_Sub_Folder_Name 'IE_Sub Record.txt'],[QQQ IE_Sub],'delimiter','\t','newline','pc','-append');
                imwrite(Test_Image,[Processed_Sub_Folder_Name sprintf('NN_%03d_%s.png',QQQ,System)],'png');
            end
        else
            if If_FolderNaming ==0
                imwrite(Lateral_Binned_Image,[Processed_Sub_Folder_Name sprintf('%03d_%s.png',QQQ,System)],'png');
            else
                imwrite(Lateral_Binned_Image,[Processed_Sub_Folder_Name folder_list(QQQ).name '.png'],'png');
            end
            if If_Save_Raw == 1
                fid = fopen([Processed_Sub_Folder_Name sprintf('%03d_%s.raw',QQQ,System)], 'w+');
                fwrite(fid, Lateral_Binned_Image', 'single');
                fclose(fid);
            end
            if If_Calculate_NN ==1
                dlmwrite([Processed_Sub_Folder_Name 'Inten_NN Record.txt'],[QQQ Intensity NN],'delimiter','\t','newline','pc','-append');
                dlmwrite([Processed_Sub_Folder_Name 'IE_Sub Record.txt'],[QQQ IE_Sub],'delimiter','\t','newline','pc','-append');
                imwrite(Test_Image,[Processed_Sub_Folder_Name sprintf('NN_%03d_%s.png',QQQ,System)],'png');
            end
        end
    else
        Processed_Sub_Folder_Name=[Data_Save_Folder sprintf('System_%s_AAF_%d_LAF_%d_RBF_%d_CBF_%d_ABTH_%g_AVE_%d_Phase_%d_ProductAVE_%d_Colomn_%d_XWidth%d',System,Axial_ave_Factor,Lateral_Ave_Factor,Row_Binning_Factor,Column_Binning_Factor,Auto_Brightness_TH,ave_factor,Phase_Shift,Product_of_Axial_Decimation_Factor_and_Ave_Factor,Colomn,X_ROI(2)-X_ROI(1)) '\'];
        if exist(Processed_Sub_Folder_Name)==0
            mkdir(Processed_Sub_Folder_Name);
        end         
        if length(Number)>0
            if If_FolderNaming ==0
                imwrite(Lateral_Binned_Image,[Processed_Sub_Folder_Name sprintf('%03d_%s.png',Number(QQQ),System)],'png');
            else
                imwrite(Lateral_Binned_Image,[Processed_Sub_Folder_Name folder_list(QQQ).name '.png'],'png');
            end
            if If_Save_Raw == 1
                fid = fopen([Processed_Sub_Folder_Name sprintf('%03d_%s.raw',QQQ,System)], 'w+');
                fwrite(fid, Lateral_Binned_Image', 'single');
                fclose(fid);
            end
            if If_Save_Raw == 1
                fid = fopen([Processed_Sub_Folder_Name sprintf('%03d_%s.raw',Number(QQQ),System)], 'w+');
                fwrite(fid, Lateral_Binned_Image', 'single');
                fclose(fid);
            end
            if If_Calculate_NN ==1
                dlmwrite([Processed_Sub_Folder_Name 'Inten_NN Record.txt'],[Number(QQQ) Intensity NN],'delimiter','\t','newline','pc','-append');
                dlmwrite([Processed_Sub_Folder_Name 'IE_Sub Record.txt'],[Number(QQQ) IE_Sub],'delimiter','\t','newline','pc','-append');
                imwrite(Test_Image,[Processed_Sub_Folder_Name sprintf('NN_%03d_%s.png',Number(QQQ),System)],'png');
            end
        elseif Y_Offset~=0
            if If_FolderNaming ==0
                imwrite(Lateral_Binned_Image,[Processed_Sub_Folder_Name sprintf('%03d_%s.png',QQQ,System)],'png');
            else
                imwrite(Lateral_Binned_Image,[Processed_Sub_Folder_Name folder_list(QQQ).name '.png'],'png');
            end
            if If_Save_Raw == 1
                fid = fopen([Processed_Sub_Folder_Name sprintf('%03d_%s.raw',QQQ,System)], 'w+');
                fwrite(fid, Lateral_Binned_Image', 'single');
                fclose(fid);
            end
            if If_Calculate_NN ==1
                dlmwrite([Processed_Sub_Folder_Name 'Inten_NN Record.txt'],[QQQ Intensity NN],'delimiter','\t','newline','pc','-append');
                dlmwrite([Processed_Sub_Folder_Name 'IE_Sub Record.txt'],[QQQ IE_Sub],'delimiter','\t','newline','pc','-append');
                imwrite(Test_Image,[Processed_Sub_Folder_Name sprintf('NN_%03d_%s.png',QQQ,System)],'png');
            end
        else
            if If_FolderNaming ==0
                imwrite(Lateral_Binned_Image,[Processed_Sub_Folder_Name sprintf('%03d_%s.png',QQQ,System)],'png');
            else
                imwrite(Lateral_Binned_Image,[Processed_Sub_Folder_Name folder_list(QQQ).name '.png'],'png');
            end
            if If_Save_Raw == 1
                fid = fopen([Processed_Sub_Folder_Name sprintf('%03d_%s.raw',QQQ,System)], 'w+');
                fwrite(fid, Lateral_Binned_Image', 'single');
                fclose(fid);
            end
            if If_Calculate_NN ==1
                dlmwrite([Processed_Sub_Folder_Name 'Inten_NN Record.txt'],[QQQ Intensity NN],'delimiter','\t','newline','pc','-append');
                dlmwrite([Processed_Sub_Folder_Name 'IE_Sub Record.txt'],[QQQ IE_Sub],'delimiter','\t','newline','pc','-append');
                imwrite(Test_Image,[Processed_Sub_Folder_Name sprintf('NN_%03d_%s.png',QQQ,System)],'png');
            end
        end
    end
%     Eff_Coef_Glass(QQQ)=max(mean(Reduced_Image(ROI_Depth_Glass(1):ROI_Depth_Glass(2),ROI_Width_Glass(1):ROI_Width_Glass(2))))
%     Eff_Coef_Sig(QQQ)=max(mean(Reduced_Image(ROI_Depth_Sig(1):ROI_Depth_Sig(2),ROI_Width_Sig(1):ROI_Width_Sig(2))))
%     Eff_Coef_BND(QQQ)=max(mean(Reduced_Image(ROI_Depth_BND(1):ROI_Depth_BND(2),ROI_Width_BND(1):ROI_Width_BND(2))))
%     
%     Eff_Coef_Glass_BND_Sub(QQQ)=max(mean(((Reduced_Image(ROI_Depth_Glass(1):ROI_Depth_Glass(2),ROI_Width_Glass(1):ROI_Width_Glass(2))).^2-Eff_Coef_BND(QQQ)^2).^0.5));
%     Eff_Coef_Sig_BND_Sub(QQQ)=max(mean(((Reduced_Image(ROI_Depth_Sig(1):ROI_Depth_Sig(2),ROI_Width_Sig(1):ROI_Width_Sig(2))).^2-Eff_Coef_BND(QQQ)).^0.5));

% %% For Noise Analysis
%     ROI_Width=[1 8];
%     ROI_Height=[1 8];
%     Noise_ROI=Reduced_Image(ROI_Height(1):ROI_Height(2),ROI_Width(1):ROI_Width(2));
%     Noise_EffMap=mean(Noise_ROI(:))
%     %imagesc(Reduced_Image);
% 
%     %colormap(gray);
    if If_Calculate_NN == 1
        Matrix_Record=[Matrix_Record;[QQQ,Intensity,NN]];
    end
 end
% subplot(3,1,1)clearvars 
% plot(ave_factor,(Eff_Coef_Glass.^2-Eff_Coef_BND.^2).^0.5);
% plot(ave_factor,Eff_Coef_Glass_BND_Sub);
% ylim([0 0.02]);
% 
% subplot(3,1,2)
% plot(ave_factor,(Eff_Coef_Sig.^2-Eff_Coef_BND.^2).^0.5);
% ylim([0 0.01]);
% 
% subplot(3,1,3)
% plot(ave_factor,Eff_Coef_Sig-Eff_Coef_BND);
% %plot(ave_factor,Eff_Coef_BND);
% ylim([0 0.01]);
if If_Calculate_NN == 1
    Test_Image=squeeze(After_Npoint_Image_Stack(X_ROI(1):X_ROI(2),NN_Y_Range,NN_Range(1):NN_Range(2)))';
    %Test_Image=squeeze(After_Npoint_Image_Stack(X_ROI(1):X_ROI(2),NN_Y_Range,:))';
    
    subplot(2,1,1)
    scatter(Intensity_Map_Noise(:),Noise_Map(:));
    xlim([0 4000]);
    ylim([0 max(Noise_Map(:))]);
    subplot(2,1,2)
    imagesc(Test_Image);
    caxis([0 10]);
    Test_Image(isnan(Test_Image))=0;
    mean_test=mean(mean(Test_Image))
end

fclose('all');

folder_list(QQQ).name
if If_Calculate_NN ==1
    NN
    hist(IE_Map,100)
end