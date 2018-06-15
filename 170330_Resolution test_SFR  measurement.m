clear all
%%
root_folder_path='D:\170329_Obj test\TIFF';



%root_folder_path='I:\161223_e2v EM1\CeYAG\Different Exp Time';
If_Best_Only=1;
Best_Index=1;
M=10;
sub_folder_path='13';

%
Camera_Pixel_Size=6.45; %micron
Adaptor_M=0.63;

Sampling_Resolution=Camera_Pixel_Size/M/Adaptor_M;

sub_processed_path='Processed Data';

Peak_Corr=1;

%Down_Sample_Ratio=50;       %Note: SFR??downsample, ???ROI
SFR_Window_Half_Width=100;


folder_path=[root_folder_path '\' sub_folder_path '\'];


Data_Name=sub_folder_path;
Display_Gain=4;


file_list=dir(folder_path);
file_list=file_list(3:end);

%


if If_Best_Only == 1
    TotalNumber=1;
    Resolution_Mean=zeros(1,1);
    Resolution_STD=zeros(1,1);
else
    TotalNumber=size(file_list,1);
    Resolution_Mean=zeros(size(file_list,1),1);
    Resolution_STD=zeros(size(file_list,1),1);

end

for Number=1:TotalNumber
    Edge_Window_Half_Size=round(50/Sampling_Resolution);
    if If_Best_Only == 1
        file_name=file_list(Best_Index).name;
    else
        file_name=file_list(Number).name;        
    end
    file_path=[folder_path file_name];

    Data_Save_Folder=[root_folder_path '\' sub_processed_path '\'];
    mkdir(Data_Save_Folder);
    Processed_Data_Path=[Data_Save_Folder Data_Name '_' file_name '.bin'];
    cd(folder_path);

    Data_Image_Raw=mean(double(imread(file_path,'tiff')),3);%.*FF;
    % 2-bound edge detection
    Left_Bound=100;
    Right_Bound=100;

    X=1:size(Data_Image_Raw,1);

    L=mean(Data_Image_Raw(:,1:Left_Bound),2);
    R=mean(Data_Image_Raw(:,(size(Data_Image_Raw,2)-Right_Bound+1):size(Data_Image_Raw,2)),2);

    L_grad=[diff(L);0];
    R_grad=[diff(R);0];
    plot(X,L,X,R,X,L_grad,X,R_grad);
    [L_Extreme_value L_Extreme]=findpeaks(abs(L_grad),'MINPEAKDISTANCE',round(80/Sampling_Resolution),'MINPEAKHEIGHT',max(abs(L_grad)/2));
    [R_Extreme_value R_Extreme]=findpeaks(abs(R_grad),'MINPEAKDISTANCE',round(80/Sampling_Resolution),'MINPEAKHEIGHT',max(abs(R_grad)/2));
    L_Extreme=L_Extreme+Peak_Corr;
    R_Extreme=R_Extreme+Peak_Corr;

    %% ?SFR target???, ???ESF????????edge??, ????????? (?1st L?2nd L)
    
    L2=mean(Data_Image_Raw(:,(Left_Bound+1):2*Left_Bound),2);
    L2_grad=[diff(L2);0];
    [L2_Extreme_value L2_Extreme]=findpeaks(abs(L2_grad),'MINPEAKDISTANCE',round(80/Sampling_Resolution),'MINPEAKHEIGHT',max(abs(L_grad)/2));

    Times=0;
    while Times<2
        if length(L_Extreme)~=length(L2_Extreme)
            if length(L_Extreme)>length(L2_Extreme)
                if abs(L_Extreme(1)-L2_Extreme(1))>50   %pixels, ?????????
                    L_Extreme(1)=[];
                elseif abs(L_Extreme(end)-L2_Extreme(end))>50
                    L_Extreme(end)=[];
                end
            elseif length(L2_Extreme)>length(L_Extreme)
                if abs(L2_Extreme(1)-L_Extreme(1))>50   %pixels
                    L2_Extreme(1)=[];
                elseif abs(L2_Extreme(end)-L_Extreme(end))>50
                    L2_Extreme(end)=[];
                end            
            end
        end
        Times=Times+1;
    end
    
    D_L_L2=Left_Bound;
    D_L_R=abs((1+Left_Bound)/2-(2*size(Data_Image_Raw,2)-Right_Bound+1)/2);
    LR_Extreme=L_Extreme+(L2_Extreme-L_Extreme)*D_L_R/D_L_L2;
    Times=0;
    
    while Times<2
        if length(L_Extreme)~=length(R_Extreme)
            if length(L_Extreme)>length(R_Extreme)
                if abs(LR_Extreme(1)-R_Extreme(1))>50   %pixels, ?????????
                    L_Extreme(1)=[];
                elseif abs(LR_Extreme(end)-R_Extreme(end))>50
                    L_Extreme(end)=[];
                end
            elseif length(R_Extreme)>length(L_Extreme)
                if abs(R_Extreme(1)-LR_Extreme(1))>50   %pixels
                    R_Extreme(1)=[];
                elseif abs(R_Extreme(end)-LR_Extreme(end))>50
                    R_Extreme(end)=[];
                end            
            end
        end
        Times=Times+1;
    end
    
    Number_of_Edge=length(L_Extreme);

    if length(R_Extreme)~=Number_of_Edge
        error('#of L&R Edge not Equal');
    end

    %% Slope calculation (for SFR)
     Angle_Array=atan(abs(L_Extreme-R_Extreme)/D_L_R)*180/pi;
     Angle_Mean=mean(Angle_Array); %degree
     Angle_STD=std(Angle_Array);
    
     

    
    %% 
    
    % Initial Guess of Edge Position
    for p=1:Number_of_Edge
        Y_Index_Array(p,1:size(Data_Image_Raw,2))=interp1([(1+Left_Bound)/2 (2*size(Data_Image_Raw,2)-Right_Bound+1)/2],[L_Extreme(p) R_Extreme(p)],[1:size(Data_Image_Raw,2)],'linear','extrap');
    end

    % Remove edges too close to boundary
    p=1;
    while(p<(Number_of_Edge+1))
        if ((max(Y_Index_Array(p,:))+Edge_Window_Half_Size)>size(Data_Image_Raw,1)) ||((min(Y_Index_Array(p,:))-Edge_Window_Half_Size)<1)
            Y_Index_Array(p,:)=[];
           p=p-1;
           Number_of_Edge=Number_of_Edge-1;
        end
        p=p+1;
    end
    Number_of_Edge=size(Y_Index_Array,1);


    % Plot Edge on Image (real)
    
    imagesc(Data_Image_Raw);
    colormap(gray);
    hold on
%     for p=1:Number_of_Edge
%         plot([(1+Left_Bound)/2,(2*size(Data_Image_Raw,2)-Right_Bound+1)/2],[L_Extreme(p),R_Extreme(p)],'Color','r','LineWidth',2)
%     end
    plot([1:size(Y_Index_Array,2)],Y_Index_Array,'.');
    hold off
    % saveas(gcf,[root_folder_path '\' sub_folder_path '.png'])

    %% Edge Window Defination (3D array for SFR, and DO NOT shift for different lateral position)

    Edge_Array=zeros(2*Edge_Window_Half_Size+1,SFR_Window_Half_Width*2+1,Number_of_Edge);
    
    for p=1:Number_of_Edge
        for q=1:(SFR_Window_Half_Width*2+1)
                Edge_Array(:,q,p)=Data_Image_Raw((round(Y_Index_Array(p,size(Data_Image_Raw,2)/2))-Edge_Window_Half_Size):(round(Y_Index_Array(p,size(Data_Image_Raw,2)/2))+Edge_Window_Half_Size),size(Data_Image_Raw,2)/2-SFR_Window_Half_Width+q);
        end
    end

    plot(Edge_Array(:,:,1)) %????????, ????pixel phase???(note??shift?, ?phase?>1 pixel)
    %% Generate pixel phase array
    Pixel_Phase_Per_Lateral_Pixel=1/tan((90-Angle_Mean)/180*pi);    %Note: unit: pixel)
    Upsampling_Ratio=floor(1/Pixel_Phase_Per_Lateral_Pixel);
    Pixel_Phase_Array=repmat((0:Pixel_Phase_Per_Lateral_Pixel:Pixel_Phase_Per_Lateral_Pixel*((SFR_Window_Half_Width*2+1)-1))',[1 Number_of_Edge]);
    %% ??????, ?????interpolation
    Manual_Adj_Phase_Ratio=1;   %????1???edge??
    Based_Index_Array_Old=[1:(2*Edge_Window_Half_Size+1)];
    Based_Index_Array_New=[1:(2*Edge_Window_Half_Size+1)*Upsampling_Ratio]/Upsampling_Ratio;
    Edge_Array_Interp=zeros(length(Based_Index_Array_New),size(Edge_Array,2),size(Edge_Array,3));
    for p=1:Number_of_Edge
        for q=1:(SFR_Window_Half_Width*2+1)
                Edge_Array_Interp(:,q,p)=interp1(Based_Index_Array_Old-Manual_Adj_Phase_Ratio*Pixel_Phase_Array(q,p),Edge_Array(:,q,p),Based_Index_Array_New); %????????, ????????
        end
    end
    %
    plot(Edge_Array_Interp(:,:,1)); %?????extrap, ?NaN????
    Edge_Array_Interp(isnan(Edge_Array_Interp))=0;
    %% Mean and Smooth
    
    Edge_Array_Interp_Mean=squeeze(mean(Edge_Array_Interp,2));
    
    %% Edge Array Normalization (Normalization AFTER SFR)
    New_Sampling_Resolution=Sampling_Resolution/Upsampling_Ratio;
    Edge_Array_Normalized=Edge_Array_Interp_Mean;
    for p=1:size(Edge_Array_Interp_Mean,2)
        Edge_Array_Normalized(:,p)=(Edge_Array_Interp_Mean(:,p)-min(Edge_Array_Interp_Mean(:,p)))/(max(Edge_Array_Interp_Mean(:,p))-min(Edge_Array_Interp_Mean(:,p)));
    end
    New_Position_Array=New_Sampling_Resolution:New_Sampling_Resolution:New_Sampling_Resolution*size(Edge_Array_Normalized,1);

    plot(New_Position_Array,Edge_Array_Normalized)

    %%
    
    
    %% ESF Resolution Calc
    Resolution_Array=zeros(1,size(Edge_Array_Normalized,2));
    for p=1:size(Edge_Array_Normalized,2)
        if Edge_Array_Normalized(1,p)>Edge_Array_Normalized(end,p)
            Resolution_Array(p)=(find(Edge_Array_Normalized(:,p)<0.1,1,'first')-find(Edge_Array_Normalized(:,p)<0.9,1,'first'))*New_Sampling_Resolution;
        else
            Resolution_Array(p)=(find(Edge_Array_Normalized(:,p)>0.9,1,'first')-find(Edge_Array_Normalized(:,p)>0.1,1,'first'))*New_Sampling_Resolution;
        end
    end
    %%
    % for p=1:size(Edge_Array_Normalized,2)
    %     plot(Edge_Array_Normalized(:,p));
    %     hold on
    % end
    plot(New_Position_Array,Edge_Array_Normalized);
    xlim([New_Position_Array(1) New_Position_Array(end)]);
    saveas(gcf,[root_folder_path '\' sub_processed_path '\' sub_folder_path sprintf('_Fig%d',Number) '_Butterfly.png'])
    Resolution_Mean(Number)=mean(Resolution_Array);
    Resolution_STD(Number)=std(Resolution_Array);
end

plot(Resolution_Mean)
xlabel('Number');
ylabel('Resolution (micron)');
saveas(gcf,[root_folder_path '\' sub_processed_path '\Resolution Chart_' sub_folder_path sprintf('_Min_%g.png',min(Resolution_Mean))])

%% Find The Best
Smooth_Window=10;
Resolution_Mean_Smooth=smooth(Resolution_Mean,Smooth_Window);

plot(Resolution_Mean_Smooth)
xlabel('Number');
ylabel('Resolution (micron)');

[Value_of_the_best_resolution Index_of_the_best_resolution]=min(Resolution_Mean_Smooth);

%% Serious calc at the best
if If_Best_Only ==1
    
end


%%
 
 NN=14;
 
 plot(Edge_Array_Normalized(:,NN));
