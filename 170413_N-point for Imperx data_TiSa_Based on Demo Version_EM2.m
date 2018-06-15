clear all
%%
Save_File_Name='170413_Imperx on EM2';
root_folder_path='D:\170413_Imperx on EM2 test\';
last_folder_name='1';
folder_path=[root_folder_path last_folder_name '\'];

SR_Lateral=0.222;
SR_Axial=0.279;

If_EffMap=0;
If_Norm=1;
Data_Save_Folder='D:\170413_Imperx on EM2 test\Processed Data\';


cd(folder_path);

% Data format related
Row=648;
Colomn=4;
Byte_Skip=0;
% Processing related
Column_Binning_Factor=1;
Y_Ave_Factor=4;
ave_factor=[16];
Product_of_Axial_Decimation_Factor_and_Ave_Factor=16;

for QQQ=1:length(ave_factor)
    N=4;
    Offset_1=0;
    Offset_2=0; 

    Max_Number_of_Frame=3000000;
    if If_EffMap ==1
        Processed_Data_Path=[Data_Save_Folder Save_File_Name last_folder_name sprintf('_EffMap_Ave_Factor_%d_CBF_%d.raw',ave_factor(QQQ),Column_Binning_Factor)];
    elseif If_EffMap ==0
        Processed_Data_Path=[Data_Save_Folder Save_File_Name last_folder_name sprintf('_Ave_Factor_%d_CBF_%d.raw',ave_factor(QQQ),Column_Binning_Factor)];
    elseif If_EffMap ==-1
        Processed_Data_Path=[Data_Save_Folder Save_File_Name last_folder_name sprintf('Test_Ave_Factor_%d_CBF_%d.raw',ave_factor(QQQ),Column_Binning_Factor)];
    end
    %%
    file_list=dir(folder_path);
    file_list=downsample(file_list(4:end),floor(Product_of_Axial_Decimation_Factor_and_Ave_Factor/ave_factor(QQQ)));
    Frame=length(file_list);
    %%
    Ave_Temp=zeros(Row,Colomn,ave_factor(QQQ));
    After_Npoint_Frame_Length=floor(Frame/N/ave_factor(QQQ));
    After_Npoint_Image_Stack=zeros(Row,Colomn/Column_Binning_Factor,min(Max_Number_of_Frame,After_Npoint_Frame_Length));
    Npoint_Temp=zeros(Row,Colomn/Column_Binning_Factor,N);

    X=[1:Frame];
    Frame_Ave=zeros([min(Max_Number_of_Frame,After_Npoint_Frame_Length)*N 1]);
    for p=1:min(Max_Number_of_Frame,After_Npoint_Frame_Length)
        for q=1:N
            for r=1:ave_factor(QQQ)

                file_path=[folder_path file_list((p-1)*N*ave_factor(QQQ)+(q-1)*ave_factor(QQQ)+r).name];

                fin=fopen(file_path);

                fseek(fin, Byte_Skip, 'bof');

                Ave_Temp(:,:,r)=fread(fin,[Row,Colomn],'uint16'); %*Frame   不知為何, 看起來就像是要除16
                fclose(fin);
            end
            Mean=mean(Ave_Temp,3);
            
            Reduced_Length=size(Mean,2)/Column_Binning_Factor;
            for tt=1:Column_Binning_Factor
               Npoint_Temp(:,:,q)=Npoint_Temp(:,:,q)+Mean(:,(Column_Binning_Factor-(tt-1)):Column_Binning_Factor:(Column_Binning_Factor*Reduced_Length)-(tt-1));
            end
            Npoint_Temp(:,:,q)=Npoint_Temp(:,:,q)/Column_Binning_Factor;
            
            
            if If_Norm ==1
                Frame_Ave((p-1)*N+q)=mean(mean(Npoint_Temp(:,:,q)));
                Npoint_Temp(:,:,q)=Npoint_Temp(:,:,q)/Frame_Ave((p-1)*N+q);
            end
        end
        
        if If_EffMap ==1
            After_Npoint_Image_Stack(:,:,p)=((N*sum(Npoint_Temp.^2,3)-sum(Npoint_Temp,3).^2).^0.5)*(2^0.5)/N./mean(Npoint_Temp,3);
        elseif If_EffMap ==0
            After_Npoint_Image_Stack(:,:,p)=((N*sum(Npoint_Temp.^2,3)-sum(Npoint_Temp,3).^2).^0.5)*(2^0.5)/N;
        elseif If_EffMap ==-1
            After_Npoint_Image_Stack(:,:,p)=mean(Npoint_Temp,3);
        end
        
        disp(p);
    end

    if If_Norm ==1
        After_Npoint_Image_Stack=After_Npoint_Image_Stack*mean(Frame_Ave);
    end
    %%
    %fid = fopen(Processed_Data_Path, 'w+');
    %fwrite(fid, After_Npoint_Image_Stack, 'double');
    %fclose(fid);
    disp(QQQ);
    
    %%
    Axial_ave_Factor=1;
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
    Reduced_Image=squeeze(mean(Reduced_Stack(:,1:Y_Ave_Factor,:),2))';
    
    Reduced_Image=Reduced_Image(1:min(size(Reduced_Image,1),Maximum_Axial_Frame),:);
    
    
%%
    clear After_Npoint_Image_Stack
    
%%

    [Max_Map Max_Index_Map]=max(Reduced_Stack,[],3);
    Mean_Map=mean(Reduced_Stack,3);
    
    %Mean_Map=mean(Reduced_Stack(:,:,(Max_Index_Map-500):(Max_Index_Map+500)),3);

    IE_Map=(Max_Map-Mean_Map)./Mean_Map;
    imagesc(IE_Map);
    
    
    %%
    
    if If_EffMap ==1
        C_max=0.01;
        C_min=0.001;
    elseif If_EffMap ==0
        if ave_factor==16
            C_max=64;
            C_min=3.2;
        elseif ave_factor==8
            C_max=32-3.2+6;
            C_min=6;
        end
    elseif If_EffMap ==-1
        C_max=32000;
        C_min=320;
    end
    Reduced_Image_normalized=(Reduced_Image-C_min)/(C_max-C_min);
    Reduced_Image_normalized(Reduced_Image_normalized<0)=0;
    Reduced_Image_normalized(Reduced_Image_normalized>1)=1;
    Reduced_Image_normalized=Reduced_Image_normalized(size(Reduced_Image_normalized,1):-1:1,:);
    subplot(1,1,1)
    imagesc(Reduced_Image_normalized);
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
    
    if If_EffMap ==1
        imwrite(Reduced_Image_normalized,[Processed_Data_Path '_Bscan_EffMap.png'],'png');
    elseif If_EffMap ==0
        imwrite(Reduced_Image_normalized,[Processed_Data_Path '_Bscan.png'],'png');
    elseif If_EffMap ==-1
        imwrite(Reduced_Image_normalized,[Processed_Data_Path '_Bscan_Test.png'],'png');
     end
%     Eff_Coef_Glass(QQQ)=max(mean(Reduced_Image(ROI_Depth_Glass(1):ROI_Depth_Glass(2),ROI_Width_Glass(1):ROI_Width_Glass(2))))
%     Eff_Coef_Sig(QQQ)=max(mean(Reduced_Image(ROI_Depth_Sig(1):ROI_Depth_Sig(2),ROI_Width_Sig(1):ROI_Width_Sig(2))))
%     Eff_Coef_BND(QQQ)=max(mean(Reduced_Image(ROI_Depth_BND(1):ROI_Depth_BND(2),ROI_Width_BND(1):ROI_Width_BND(2))))
%     
%     Eff_Coef_Glass_BND_Sub(QQQ)=max(mean(((Reduced_Image(ROI_Depth_Glass(1):ROI_Depth_Glass(2),ROI_Width_Glass(1):ROI_Width_Glass(2))).^2-Eff_Coef_BND(QQQ)^2).^0.5));
%     Eff_Coef_Sig_BND_Sub(QQQ)=max(mean(((Reduced_Image(ROI_Depth_Sig(1):ROI_Depth_Sig(2),ROI_Width_Sig(1):ROI_Width_Sig(2))).^2-Eff_Coef_BND(QQQ)).^0.5));

%% For Noise Analysis
    ROI_Width=[1 8];
    ROI_Height=[1 8];
    Noise_ROI=Reduced_Image(ROI_Height(1):ROI_Height(2),ROI_Width(1):ROI_Width(2));
    Noise_EffMap=mean(Noise_ROI(:))
    %imagesc(Reduced_Image);

    %colormap(gray);
end
% subplot(3,1,1)
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

fclose('all');

