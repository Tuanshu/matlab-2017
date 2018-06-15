clearvars 
%%

Save_File_Name='170427_Hilbert on FF-OCT';
root_folder_path='C:\TuanShu\170427_Hilbert on FF-OCT\';
last_folder_name=[26];


If_EM2=1;
If_Auto_Brightness=1;
Auto_Brightness_TH=0.99;
N=4;

If_EffMap=0;%0;
If_Norm=1;

If_Enhace_Deep_Signal=1;


Data_Save_Folder=[root_folder_path '\Processed Data\'];

If_AlterAxialProcess=1;

% Data format related
Row=648;
Colomn=8;
Page=319527;
Byte_Skip=0;
% Processing related
Column_Binning_Factor=1;
Row_Binning_Factor=2;
Lateral_Ave_Factor=1;   %New Param 170413
Axial_Ave_Factor=8*4; %8*8=64, *4 for 4-point

%Y_Ave_Factor=8;
Product_of_Axial_Decimation_Factor_and_Ave_Factor=round(8);    %should be 64 for 0.28/4, use 8 here
ave_factor=Product_of_Axial_Decimation_Factor_and_Ave_Factor;%[4*4*4];

Phase_Shift=0;

for QQQ=1:length(last_folder_name)
    folder_path=[root_folder_path num2str(last_folder_name(QQQ)) '\'];
    cd(folder_path);
    micron_per_frame=0.2/ave_factor/N;
    Offset_1=0;
    Offset_2=0; 

    Max_Number_of_Frame=3000000;
    if If_EffMap ==1
        Processed_Data_Path=[Data_Save_Folder Save_File_Name num2str(last_folder_name(QQQ))];
    elseif If_EffMap ==0
        Processed_Data_Path=[Data_Save_Folder Save_File_Name num2str(last_folder_name(QQQ))];
    elseif If_EffMap ==-1
        Processed_Data_Path=[Data_Save_Folder Save_File_Name num2str(last_folder_name(QQQ))];
    end
    %%
    file_list=dir(folder_path);    
    for k = length(file_list):-1:1
        if file_list(k).isdir
            file_list(k) = [ ];
            continue
        end
        fname = file_list(k).name;
        if (fname(1) == '.') || strcmp(fname,'Processed Data') ||  strcmp(fname,'Other')
            file_list(k) = [ ];
        end
    end
    %%
    if length(file_list)==1
        fin = fopen(file_list(1).name);
        Image_Temp=fread(fin,[Row,inf],'uint16','b');
        fclose(fin);
        
        %
        Colomn_Total=size(Image_Temp,2);

        Frame=Colomn_Total/Colomn;
        Image_Stack=zeros(Row,Colomn,Frame);
        for r=1:Frame
            Image_Stack(:,:,r)=Image_Temp(:,(1+(r-1)*Colomn):(r*Colomn));
        end
        
        clear Image_Temp
        
    else
        error('Multiple Files in Folder');
    end
    
    %%
    Frame_Ave=(mean(mean(Image_Stack,1),2));
    ADJ_COEF=mean(Frame_Ave)./Frame_Ave;
    if If_Norm
        Image_Stack=Image_Stack.*repmat(Frame_Ave,size(Image_Stack,1),size(Image_Stack,2),1);
    end
  %% Row Binning
  
    Temp=0;
    Reduced_Width=size(Image_Stack,1)/Row_Binning_Factor;
    for tt=1:Row_Binning_Factor
       Temp=Temp+Image_Stack((Row_Binning_Factor-(tt-1)):Row_Binning_Factor:(Row_Binning_Factor*Reduced_Width)-(tt-1),:,:);
    end
    Image_Stack_Row_Binned=Temp/Row_Binning_Factor;
    
    clear Image_Stack
    %% Axial Binning (replacing previous binning before N-point))
    Temp=0;
    for p=1:ave_factor
        Temp=Temp+Image_Stack_Row_Binned(:,:,(ave_factor-(p-1)):ave_factor:(ave_factor*floor(size(Image_Stack_Row_Binned,3)/ave_factor))-(p-1));
    end
    
    Image_Stack_Axial_Binned=Temp/ave_factor;
    
    clear Image_Stack Temp
    %% FFT
    FFT_Stack=fft(Image_Stack_Axial_Binned,[],3);
    
    Normalized_F=(1/size(Image_Stack_Axial_Binned,3):1/size(Image_Stack_Axial_Binned,3):1)/ave_factor;    %Normalized freq based on original frame rate
    plot(Normalized_F,squeeze(real(FFT_Stack(300,4,:))));
    
    %clear Image_Stack_Axial_Binned
%% Filtering
    SB_Norm=0.0035;
    LB_Norm=0.006;
    
    FFT_Stack_Filtered=FFT_Stack;
    FFT_Stack_Filtered(:,:,1:ceil(SB_Norm*size(FFT_Stack,3)*ave_factor))=0;
    FFT_Stack_Filtered(:,:,floor(LB_Norm*size(FFT_Stack,3)*ave_factor):end)=0;
    
    %clear FFT_Stack
    
    Image_Stack_Hilbert=2*abs(ifft(FFT_Stack_Filtered,[],3));
    
    clear FFT_Stack_Filtered

    disp(QQQ);

    %
    Maximum_Axial_Frame=20000;
    Temp=0;
    Axial_Length_Original=size(Image_Stack_Hilbert,3);
    Axial_Length_Used=floor(Axial_Length_Original/Axial_Ave_Factor)*Axial_Ave_Factor;
    Reduced_Length=Axial_Length_Used/Axial_Ave_Factor;
    for p=1:Axial_Ave_Factor
       Temp=Temp+Image_Stack_Hilbert(:,:,(Axial_Ave_Factor-(p-1)):Axial_Ave_Factor:(Axial_Ave_Factor*Reduced_Length)-(p-1));
    end
    Reduced_Stack=Temp/Axial_Ave_Factor;
    
    %Reduced_Image=squeeze(Reduced_Stack(:,size(Reduced_Stack,2)/2,:))';
    %Reduced_Image=squeeze(mean(Reduced_Stack(:,1:Y_Ave_Factor/Column_Binning_Factor,:),2))';
    Reduced_Image=squeeze(mean(Reduced_Stack(:,:,:),2))';
    disp(size(Reduced_Stack,2))
    Reduced_Image=Reduced_Image(1:min(size(Reduced_Image,1),Maximum_Axial_Frame),:);
    if If_EM2 ==1
        Reduced_Image=Reduced_Image(size(Reduced_Image,1):-1:1,:);
    end
        %
    clear Image_Stack_Hilbert
    

    % Deep Signal Analysis
    if If_Enhace_Deep_Signal
        BI=1;
        BF=0.1;
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
    % Hist for auto brightness

    if If_Auto_Brightness == 1
        Histogram=hist(Reduced_Image(:),max(Reduced_Image(:)));

        Total=sum(Histogram);
        Int_Hist = cumsum(Histogram)/Total;

        %plot(Int_Hist);

        C_max_Auto=find(Int_Hist>Auto_Brightness_TH,1,'first');
        C_max=C_max_Auto
        C_min=C_max*0.05;
    else

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
    end
    
    %
    Reduced_Image_normalized=(Reduced_Image-C_min)/(C_max-C_min);
    Reduced_Image_normalized(Reduced_Image_normalized<0)=0;
    Reduced_Image_normalized(Reduced_Image_normalized>1)=1;
    
    
    % Lateral Ave
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
    %axis equal
    caxis([0 1]);
    set(gca,'xtick',[0  20  40  60  80  100],'xticklabel',[0  20  40  60  80  100]*0.33,'ytick',[0  100  200  300  400  500 600],'yticklabel',[0  100  200  300  400  500 600]*0.2645);
    xlabel('Lateral Position (micron)');
    ylabel('Axial Position (micron)');
    colormap(gray);
    
    if If_EffMap ==1
        imwrite(Lateral_Binned_Image,[Processed_Data_Path '_Bscan_EffMap.png'],'png');
    elseif If_EffMap ==0
        if If_Enhace_Deep_Signal
            imwrite(Lateral_Binned_Image,[Processed_Data_Path sprintf('_Bscan_AAF_%d_LAF_%d_CBF_%d_ABTH_%g_AVE_%d_Phase_%d_ProductAVE_%d_Colomn_%d_Hilbert_SB%g_LB%g_DeepEnhace.png',Axial_Ave_Factor,Lateral_Ave_Factor,Column_Binning_Factor,Auto_Brightness_TH,ave_factor,Phase_Shift,Product_of_Axial_Decimation_Factor_and_Ave_Factor,Colomn,SB_Norm,LB_Norm)],'png');
        else
            imwrite(Lateral_Binned_Image,[Processed_Data_Path sprintf('_Bscan_AAF_%d_LAF_%d_CBF_%d_ABTH_%g_AVE_%d_Phase_%d_ProductAVE_%d_Colomn_%d_Hilbert_SB%g_LB%g.png',Axial_Ave_Factor,Lateral_Ave_Factor,Column_Binning_Factor,Auto_Brightness_TH,ave_factor,Phase_Shift,Product_of_Axial_Decimation_Factor_and_Ave_Factor,Colomn,SB_Norm,LB_Norm)],'png');
        end
    elseif If_EffMap ==-1
        imwrite(Lateral_Binned_Image,[Processed_Data_Path '_Bscan_Test.png'],'png');
    end
end

fclose('all');

