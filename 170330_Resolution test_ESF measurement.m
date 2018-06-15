clear all
%%
root_folder_path='D:\170329_Obj test\TIFF';



%root_folder_path='I:\161223_e2v EM1\CeYAG\Different Exp Time';
If_Best_Only=1;
Best_Index=33;
M=10;
sub_folder_path='10';

%
Camera_Pixel_Size=6.45; %micron
Adaptor_M=0.63;

Sampling_Resolution=Camera_Pixel_Size/M/Adaptor_M;

sub_processed_path='Processed Data';

Peak_Corr=1;

Down_Sample_Ratio=50;


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
    Edge_Window_Half_Size=50/Sampling_Resolution;
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

    % 
    % [L_Maxes_value L_Maxes]=findpeaks(L_grad,'MINPEAKDISTANCE',50);
    % [L_Mins_value L_Mins]=findpeaks(max(L_grad)-L_grad,'MINPEAKDISTANCE',50);
    % [R_Maxes_value R_Maxes]=findpeaks(R_grad,'MINPEAKDISTANCE',50);
    % [R_Mins_value R_Mins]=findpeaks(max(R_grad)-R_grad,'MINPEAKDISTANCE',50);

    if length(L_Extreme)~=length(R_Extreme)
        if length(L_Extreme)>length(R_Extreme)
            if abs(L_Extreme(1)-R_Extreme(1))>50   %pixels
                L_Extreme(1)=[];
            elseif abs(L_Extreme(end)-R_Extreme(end))>50
                L_Extreme(end)=[];
            end
        elseif length(R_Extreme)>length(L_Extreme)
            if abs(R_Extreme(1)-L_Extreme(1))>50   %pixels
                R_Extreme(1)=[];
            elseif abs(R_Extreme(end)-L_Extreme(end))>50
                R_Extreme(end)=[];
            end            
        end
    end
    
    if length(L_Extreme)~=length(R_Extreme)     %Do this twice
        if length(L_Extreme)>length(R_Extreme)
            if abs(L_Extreme(1)-R_Extreme(1))>50   %pixels
                L_Extreme(1)=[];
            elseif abs(L_Extreme(end)-R_Extreme(end))>50
                L_Extreme(end)=[];
            end
        elseif length(R_Extreme)>length(L_Extreme)
            if abs(R_Extreme(1)-L_Extreme(1))>50   %pixels
                R_Extreme(1)=[];
            elseif abs(R_Extreme(end)-L_Extreme(end))>50
                R_Extreme(end)=[];
            end            
        end
    end
    
    
    
    Number_of_Edge=length(L_Extreme);

    if length(R_Extreme)~=Number_of_Edge
        error('#of L&R Edge not Equal');
    end

    % Initial Guess of Edge Position
    for p=1:Number_of_Edge
        Y_Index_Array(p,1:size(Data_Image_Raw,2))=interp1([1 size(Data_Image_Raw,2)],[L_Extreme(p) R_Extreme(p)],[1:size(Data_Image_Raw,2)]);
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


    % Plot Edge on Image
    
%     imagesc(Data_Image_Raw);
%     colormap(gray);
%     hold on
%     for p=1:Number_of_Edge
%         plot([1,size(Data_Image_Raw,2)],[L_Extreme(p),R_Extreme(p)],'Color','r','LineWidth',2)
%     end
%     hold off
    % saveas(gcf,[root_folder_path '\' sub_folder_path '.png'])

    %% Edge Window Defination

    Edge_Array=zeros(2*Edge_Window_Half_Size+1,floor(size(Data_Image_Raw,2)/Down_Sample_Ratio)*Number_of_Edge);

    for p=1:Number_of_Edge
        for q=1:floor(size(Data_Image_Raw,2)/Down_Sample_Ratio)
                Edge_Array(:,(p-1)*floor(size(Data_Image_Raw,2)/Down_Sample_Ratio)+q)=Data_Image_Raw(round((Y_Index_Array(p,q)-Edge_Window_Half_Size):(Y_Index_Array(p,q)+Edge_Window_Half_Size)),(q-1)*Down_Sample_Ratio+1);
        end
    end

    %% Edge Array Normalization
    Edge_Array_Normalized=Edge_Array;
    for p=1:size(Edge_Array,2)
        Edge_Array_Normalized(:,p)=(Edge_Array(:,p)-min(Edge_Array(:,p)))/(max(Edge_Array(:,p))-min(Edge_Array(:,p)));
    end

    %% Resolution Calc
    Resolution_Array=zeros(1,size(Edge_Array,2));
    for p=1:size(Edge_Array,2)
        if Edge_Array_Normalized(1,p)>Edge_Array_Normalized(end,p)
            Resolution_Array(p)=(find(Edge_Array_Normalized(:,p)<0.1,1,'first')-find(Edge_Array_Normalized(:,p)<0.9,1,'first'))*Sampling_Resolution;
        else
            Resolution_Array(p)=(find(Edge_Array_Normalized(:,p)>0.9,1,'first')-find(Edge_Array_Normalized(:,p)>0.1,1,'first'))*Sampling_Resolution;
        end
    end
    %%
    % for p=1:size(Edge_Array_Normalized,2)
    %     plot(Edge_Array_Normalized(:,p));
    %     hold on
    % end
    plot(Edge_Array_Normalized);
    xlim([0 Edge_Window_Half_Size*2+1]);
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
