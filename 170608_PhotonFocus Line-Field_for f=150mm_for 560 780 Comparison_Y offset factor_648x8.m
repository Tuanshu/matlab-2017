clearvars 
%%

Save_File_Name='170612_EM2 780 with PhotonFocus';
root_folder_path='F:\170612_PhotonFocus Test\';
%last_folder_name=[73:100];%
%last_folder_name=[1:3 5:6 12 17:20 22:32 34:38 40:49 51:56 58:59 61:63 66 68:69 71 74:86 88:91 93 96:97 100];%
last_folder_name='128x_2080fps_Forearm_648x8';%
Number=[1];

System='780';

If_EM2=1;
If_Auto_Brightness=1;
Auto_Brightness_TH=0.99;
N=4;

If_Norm=1;
Data_Save_Folder=[root_folder_path '\Processed Data\'];
if exist(Data_Save_Folder)==0
    mkdir(Data_Save_Folder);
end
If_AlterAxialProcess=1;
If_Enhace_Deep_Signal=1;
BI=1;
BF=0.2;
% Data format related
Row=648;
Colomn=8;
Y_Offset=0;
Byte_Skip=Y_Offset*Row*2;
% Processing related
Column_Binning_Factor=1;
Row_Binning_Factor=1;
Lateral_Ave_Factor=1;   %New Param 170413

if strcmp(System,'560')
    Axial_ave_Factor=4;
elseif strcmp(System,'780')
    Axial_ave_Factor=2;
end

%Y_Ave_Factor=8;
Product_of_Axial_Decimation_Factor_and_Ave_Factor=round(1);
ave_factor=Product_of_Axial_Decimation_Factor_and_Ave_Factor;%[4*4*4];

Phase_Shift=0;

for QQQ=1:length(Number)
    folder_path=[root_folder_path sprintf('%s_%g',last_folder_name,Number(QQQ)) '\'];
    cd(folder_path);
    micron_per_frame=0.2/ave_factor/N;
    Offset_1=0;
    Offset_2=0; 

    Max_Number_of_Frame=3000000;
    Processed_Data_Path=[Data_Save_Folder Save_File_Name '_' num2str(Number(QQQ))];

    %%
    file_list=dir(folder_path);
    file_list=downsample(file_list(4:end),floor(Product_of_Axial_Decimation_Factor_and_Ave_Factor/ave_factor),Phase_Shift);
    Frame=length(file_list);
    %%
    Ave_Temp=zeros(Row,Colomn,ave_factor);
    After_Npoint_Frame_Length=floor(Frame/N/ave_factor);
    After_Npoint_Image_Stack=zeros(Row/Row_Binning_Factor,Colomn/Column_Binning_Factor,min(Max_Number_of_Frame,After_Npoint_Frame_Length));

    
    Npoint_Temp=zeros(Row/Row_Binning_Factor,Colomn/Column_Binning_Factor,N);

    X=[1:Frame];
    Frame_Ave=zeros([min(Max_Number_of_Frame,After_Npoint_Frame_Length)*N 1]);
    for p=1:min(Max_Number_of_Frame,After_Npoint_Frame_Length)
        for q=1:N
            for r=1:ave_factor
%%
                file_path=[folder_path file_list((p-1)*N*ave_factor+(q-1)*ave_factor+r).name];

                fin=fopen(file_path);

                fseek(fin, Byte_Skip, 'bof');
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
            
            
            
            if If_Norm ==1
                Frame_Ave((p-1)*N+q)=mean(mean(Npoint_Temp(:,:,q)));
                Npoint_Temp(:,:,q)=Npoint_Temp(:,:,q)/Frame_Ave((p-1)*N+q);
            end
        end
        
        After_Npoint_Image_Stack(:,:,p)=((N*sum(Npoint_Temp.^2,3)-sum(Npoint_Temp,3).^2).^0.5)*(2^0.5)/N;

        
        
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

    if If_EM2 ==1
        Reduced_Image=Reduced_Image(size(Reduced_Image,1):-1:1,:);
    end   
    
    %%
    
    if If_Enhace_Deep_Signal

        Linear_Axial_Profile=[BI:-1*(BI-BF)/(size(Reduced_Image,1)-1):BF]';
%         imagesc(Reduced_Image)
%         caxis([C_min C_max])

%         Profile_Smoothing_Window=10;
% 
%         Image_Axial_Profile=mean(Reduced_Image,2);
% 
%         Image_Axial_Profile_Smooth=conv(Image_Axial_Profile,ones(Profile_Smoothing_Window,1)/Profile_Smoothing_Window,'same');
%         plot(Image_Axial_Profile_Smooth)
        Enhace_Coef=mean(Linear_Axial_Profile)./Linear_Axial_Profile;
        
        Reduced_Image=Reduced_Image.*repmat(Enhace_Coef,1,size(Reduced_Image,2));
        
    end    
%%        
    clear After_Npoint_Image_Stack
    

     
    
    if If_Auto_Brightness == 1
        Histogram=hist(Reduced_Image(:),max(Reduced_Image(:)));

        Total=sum(Histogram);
        Int_Hist = cumsum(Histogram)/Total;

        plot(Int_Hist);

        C_max_Auto=find(Int_Hist>Auto_Brightness_TH,1,'first');
        C_max=C_max_Auto
        C_min=C_max*0.05;
    else
    %% Hist for auto brightness



    if ave_factor==16
        C_max=32;
        C_min=3.2;
    elseif ave_factor==8
        C_max=32-3.2+6;
        C_min=6;
    end
    
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
    

    if If_Enhace_Deep_Signal
        Processed_Sub_Folder_Name=[Data_Save_Folder sprintf('System_%s_AAF_%d_LAF_%d_RBF_%d_CBF_%d_ABTH_%g_AVE_%d_Phase_%d_ProductAVE_%d_Colomn_%d_DeepEnhace_BI%g_BF%g',System,Axial_ave_Factor,Lateral_Ave_Factor,Row_Binning_Factor,Column_Binning_Factor,Auto_Brightness_TH,ave_factor,Phase_Shift,Product_of_Axial_Decimation_Factor_and_Ave_Factor,Colomn,BI,BF) '\'];
        if exist(Processed_Sub_Folder_Name)==0
            mkdir(Processed_Sub_Folder_Name);
        end
        if length(Number)>1
            imwrite(Lateral_Binned_Image,[Processed_Sub_Folder_Name sprintf('%03d_%s.png',QQQ,System)],'png');
        elseif Y_Offset~=0
            imwrite(Lateral_Binned_Image,[Processed_Sub_Folder_Name sprintf('%03d_%s_Y%d.png',Number(QQQ),System,Y_Offset)],'png');
        else
            imwrite(Lateral_Binned_Image,[Processed_Sub_Folder_Name sprintf('%03d_%s.png',Number(QQQ),System)],'png');
        end
    else
        Processed_Sub_Folder_Name=[Data_Save_Folder sprintf('System_%s_AAF_%d_LAF_%d_RBF_%d_CBF_%d_ABTH_%g_AVE_%d_Phase_%d_ProductAVE_%d_Colomn_%d',System,Axial_ave_Factor,Lateral_Ave_Factor,Row_Binning_Factor,Column_Binning_Factor,Auto_Brightness_TH,ave_factor,Phase_Shift,Product_of_Axial_Decimation_Factor_and_Ave_Factor,Colomn) '\'];
        if exist(Processed_Sub_Folder_Name)==0
            mkdir(Processed_Sub_Folder_Name);
        end         
        if length(Number)>1
            imwrite(Lateral_Binned_Image,[Processed_Sub_Folder_Name sprintf('%03d_%s.png',QQQ,System)],'png');
        elseif Y_Offset~=0
            imwrite(Lateral_Binned_Image,[Processed_Sub_Folder_Name sprintf('%03d_%s_Y%d.png',Number(QQQ),System,Y_Offset)],'png');
        else
            imwrite(Lateral_Binned_Image,[Processed_Sub_Folder_Name sprintf('%03d_%s.png',Number(QQQ),System)],'png');
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

fclose('all');

