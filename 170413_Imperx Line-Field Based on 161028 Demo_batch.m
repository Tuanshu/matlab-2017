clear all
%%
Save_File_Name='170413_Review Imperx Data';
root_folder_path='F:\170413_Review 161028 Data_Imperx\';
last_folder_name=[13];

If_EM2=0;

If_EffMap=0;
If_Norm=1;
Data_Save_Folder='F:\170413_Review 161028 Data_Imperx\Processed Data\';



% Data format related
Row=648;
Colomn=4;
Byte_Skip=0;
% Processing related
Column_Binning_Factor=4;
Row_Binning_Factor=1;   %New Param 170413
Y_Ave_Factor=4;
ave_factor=[16];
Product_of_Axial_Decimation_Factor_and_Ave_Factor=16;

for QQQ=1:length(last_folder_name)
    folder_path=[root_folder_path num2str(last_folder_name(QQQ)) '\'];
    cd(folder_path);
    N=4;
    micron_per_frame=0.2/ave_factor/N;
    Offset_1=0;
    Offset_2=0; 

    Max_Number_of_Frame=3000000;
    if If_EffMap ==1
        Processed_Data_Path=[Data_Save_Folder Save_File_Name num2str(last_folder_name(QQQ)) sprintf('_EffMap_Ave_Factor_%d_CBF_%d.raw',ave_factor,Column_Binning_Factor)];
    elseif If_EffMap ==0
        Processed_Data_Path=[Data_Save_Folder Save_File_Name num2str(last_folder_name(QQQ)) sprintf('_Ave_Factor_%d_CBF_%d.raw',ave_factor,Column_Binning_Factor)];
    elseif If_EffMap ==-1
        Processed_Data_Path=[Data_Save_Folder Save_File_Name num2str(last_folder_name(QQQ)) sprintf('Test_Ave_Factor_%d_CBF_%d.raw',ave_factor,Column_Binning_Factor)];
    end
    %%
    file_list=dir(folder_path);
    file_list=downsample(file_list(4:end),floor(Product_of_Axial_Decimation_Factor_and_Ave_Factor/ave_factor));
    Frame=length(file_list);
    %%
    Ave_Temp=zeros(Row,Colomn,ave_factor);
    After_Npoint_Frame_Length=floor(Frame/N/ave_factor);
    After_Npoint_Image_Stack=zeros(Row,Colomn/Column_Binning_Factor,min(Max_Number_of_Frame,After_Npoint_Frame_Length));
    Npoint_Temp=zeros(Row,Colomn/Column_Binning_Factor,N);

    X=[1:Frame];
    Frame_Ave=zeros([min(Max_Number_of_Frame,After_Npoint_Frame_Length)*N 1]);
    for p=1:min(Max_Number_of_Frame,After_Npoint_Frame_Length)
        for q=1:N
            for r=1:ave_factor

                file_path=[folder_path file_list((p-1)*N*ave_factor+(q-1)*ave_factor+r).name];

                fin=fopen(file_path);

                fseek(fin, Byte_Skip, 'bof');

                Ave_Temp(:,:,r)=fread(fin,[Row,Colomn],'uint16'); %*Frame   不知為何, 看起來就像是要除16
                %fclose(fin);
                fclose('all');

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
    %Reduced_Image=squeeze(mean(Reduced_Stack(:,1:Y_Ave_Factor/Column_Binning_Factor,:),2))';
    Reduced_Image=squeeze(mean(Reduced_Stack(:,:,:),2))';
    disp(size(Reduced_Stack,2))
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
    
    ROI_Depth_Glass=[1 5];
    ROI_Width_Glass=[200 250];
    
    ROI_Depth_Sig=[111 115];
    ROI_Width_Sig=[200 250];
  
    ROI_Depth_BND=[21 25];
    ROI_Width_BND=[200 250];
    
    if If_EffMap ==1
        C_max=0.01;
        C_min=0.001;
    elseif If_EffMap ==0
        if ave_factor==16
            C_max=32;
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
    
    if If_EM2 ==1
        Reduced_Image_normalized=Reduced_Image_normalized(size(Reduced_Image_normalized,1):-1:1,:);
    end
    
    %% Row Binning
    Row_Binning_Factor=1;
    Temp=0;
    Lateral_Length_Original=size(Reduced_Image_normalized,2);
    Reduced_Length=Lateral_Length_Original/Row_Binning_Factor;
    for p=1:Row_Binning_Factor
       Temp=Temp+Reduced_Image_normalized(:,(Row_Binning_Factor-(p-1)):Row_Binning_Factor:(Row_Binning_Factor*Reduced_Length)-(p-1));
    end
    Lateral_Binned_Image=Temp/Row_Binning_Factor;
    
    %Reduced_Image=squeeze(Reduced_Stack(:,size(Reduced_Stack,2)/2,:))';
    %Reduced_Image=squeeze(mean(Reduced_Stack(:,1:Y_Ave_Factor,:),2))';
    
    
    
    
    %
    
    subplot(1,1,1)
    imagesc(Lateral_Binned_Image);
    axis off
    %axis equal
    caxis([0 1]);
    set(gca,'xtick',[0  20  40  60  80  100],'xticklabel',[0  20  40  60  80  100]*0.33,'ytick',[0  100  200  300  400  500 600],'yticklabel',[0  100  200  300  400  500 600]*0.2645);
    xlabel('Lateral Position (micron)');
    ylabel('Axial Position (micron)');
    colormap(gray);
    
    if If_EffMap ==1
        imwrite(Lateral_Binned_Image,[Processed_Data_Path '_Bscan_EffMap.png'],'png');
    elseif If_EffMap ==0
        imwrite(Lateral_Binned_Image,[Processed_Data_Path sprintf('_Bscan_AAF_%d_RBF_%d_CBF_%d.png',Axial_ave_Factor,Row_Binning_Factor,Column_Binning_Factor)],'png');
    elseif If_EffMap ==-1
        imwrite(Lateral_Binned_Image,[Processed_Data_Path '_Bscan_Test.png'],'png');
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

