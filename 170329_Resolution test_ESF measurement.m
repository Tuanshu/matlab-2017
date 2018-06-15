clear all
%%
root_folder_path='D:\170329_Obj test\TIFF';

%root_folder_path='I:\161223_e2v EM1\CeYAG\Different Exp Time';

sub_folder_path='37';
sub_FF_path='FF';

Peak_Corr=1;


folder_path=[root_folder_path '\' sub_folder_path '\'];

FF_path=[root_folder_path '\' sub_FF_path '\'];

Data_Name=sub_folder_path;
Display_Gain=4;


file_list=dir(folder_path);
file_list=file_list(3:end);

FF_list=dir(FF_path);
FF_list=FF_list(3:end);
% FF
FF_name=FF_list(1).name;
FF_Raw=mean(double(imread([FF_path FF_name],'tiff')),3);

FF=1./FF_Raw;
FF=FF./max(FF(:));

%

Number=4;
for Number=1:size(file_list,1)
    Edge_Window_Half_Size=50;
    file_name=file_list(Number).name;
    file_path=[folder_path file_name];

    Data_Save_Folder='D:\170328_Resolution test\Processed Data\';
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
    [L_Extreme_value L_Extreme]=findpeaks(abs(L_grad),'MINPEAKDISTANCE',50,'MINPEAKHEIGHT',max(abs(L_grad)/2));
    [R_Extreme_value R_Extreme]=findpeaks(abs(R_grad),'MINPEAKDISTANCE',50,'MINPEAKHEIGHT',max(abs(R_grad)/2));
    L_Extreme=L_Extreme+Peak_Corr;
    R_Extreme=R_Extreme+Peak_Corr;

    % 
    % [L_Maxes_value L_Maxes]=findpeaks(L_grad,'MINPEAKDISTANCE',50);
    % [L_Mins_value L_Mins]=findpeaks(max(L_grad)-L_grad,'MINPEAKDISTANCE',50);
    % [R_Maxes_value R_Maxes]=findpeaks(R_grad,'MINPEAKDISTANCE',50);
    % [R_Mins_value R_Mins]=findpeaks(max(R_grad)-R_grad,'MINPEAKDISTANCE',50);

    Number_of_Edge=length(L_Extreme);

    if length(R_Extreme)~=Number_of_Edge
        error('#of L&R Edge not Equal');
    end

    % Initial Guess of Edge Position
    for p=1:Number_of_Edge
        Y_Index_Array(p,1:size(Data_Image_Raw,2))=interp1([1 size(Data_Image_Raw,2)],[L_Extreme(p) R_Extreme(p)],[1:size(Data_Image_Raw,2)]);
    end

    % Remove edges too close to boundary
    for p=1:Number_of_Edge
        if ((max(Y_Index_Array(p,:))+Edge_Window_Half_Size)>size(Data_Image_Raw,1)) ||((min(Y_Index_Array(p,:))-Edge_Window_Half_Size)<1)
            Y_Index_Array(p,:)=[];
        end
    end
    Number_of_Edge=size(Y_Index_Array,1);


    % Plot Edge on Image
    % 
    % imagesc(Data_Image_Raw);
    % colormap(gray);
    % hold on
    % for p=1:Number_of_Edge
    %     plot([1,size(Data_Image_Raw,2)],[L_Extreme(p),R_Extreme(p)],'Color','r','LineWidth',2)
    % end
    % hold off
    % saveas(gcf,[root_folder_path '\' sub_folder_path '.png'])

    %% Edge Window Defination

    Edge_Array=zeros(2*Edge_Window_Half_Size+1,size(Data_Image_Raw,2)*Number_of_Edge);

    for p=1:Number_of_Edge
        for q=1:size(Data_Image_Raw,2)
                Edge_Array(:,(p-1)*size(Data_Image_Raw,2)+q)=Data_Image_Raw(round((Y_Index_Array(p,q)-Edge_Window_Half_Size):(Y_Index_Array(p,q)+Edge_Window_Half_Size)),q);
        end
    end

    %% Edge Array Normalization
    Edge_Array_Normalized=Edge_Array;
    for p=1:size(Edge_Array,2)
        Edge_Array_Normalized(:,p)=(Edge_Array(:,p)-min(Edge_Array(:,p)))/(max(Edge_Array(:,p))-min(Edge_Array(:,p)));
    end

    %%
    % for p=1:size(Edge_Array_Normalized,2)
    %     plot(Edge_Array_Normalized(:,p));
    %     hold on
    % end
    plot(Edge_Array_Normalized);
    xlim([0 Edge_Window_Half_Size*2+1]);
    saveas(gcf,[root_folder_path '\' sub_folder_path sprintf('_Fig%d',Number) 'Butterfly.png'])
end
%% Supremum
% Supremum_Filter_Size_Half_Size=150;
% Supremum_Filter=fspecial('gaussian', Supremum_Filter_Size_Half_Size*2+1,(Supremum_Filter_Size_Half_Size*2+1)/3);    %gaussian mask
% Data_Image_Sub=zeros(size(Data_Image_Raw,1),size(Data_Image_Raw,2));
% for p=1:size(Data_Image_Raw,1)
%     for q=1:size(Data_Image_Raw,2)
%         x_min=max((p-Supremum_Filter_Size_Half_Size+1),1);
%         x_max=min((p+Supremum_Filter_Size_Half_Size),size(Data_Image_Raw,1));
%         y_min=max((q-Supremum_Filter_Size_Half_Size+1),1);
%         y_max=min((q+Supremum_Filter_Size_Half_Size),size(Data_Image_Raw,2));
%         Window_Temp=Supremum_Filter((x_min-p+Supremum_Filter_Size_Half_Size+1):(x_max-(p+Supremum_Filter_Size_Half_Size)+Supremum_Filter_Size_Half_Size*2+1),(y_min-q+Supremum_Filter_Size_Half_Size+1):(y_max-(q+Supremum_Filter_Size_Half_Size)+Supremum_Filter_Size_Half_Size*2+1)).*Data_Image_Raw(x_min:x_max,y_min:y_max);
%         Data_Image_Sub(p,q)=max(Window_Temp(:));
%     end
%     disp(p)
% end

% %%
% Ball_Size=10;
% Element=fspecial('disk', Ball_Size*3);
% Element(Element>0)=1;
% Image_Open=imopen(Data_Image_Raw,Element);
% Image_Close=imclose(Data_Image_Raw,Element);
% subplot(2,1,1)
% imagesc(Image_Open);
% subplot(2,1,2)
% imagesc(Image_Close);
% %% Grad
% Ball_Size=50;
% %Ball=fspecial('gaussian', Ball_Size,Ball_Size/3);    %gaussian mask
% Sobel=fspecial('Sobel');    %gaussian mask
% 
% %Data_Image_Smooth=(filter2(Ball,Data_Image_Raw,'same'));
% Data_Image_Edge_Enhaced=(filter2(Sobel,Data_Image_Raw,'same'));
% 
% imagesc(Data_Image_Edge_Enhaced);
% Data_Image_Smooth=(filter2(Ball,Data_Image_Edge_Enhaced,'same'));
% 
% % Edge Detection
% TH=0.01;
% 
% Edge_Map = edge(Data_Image_Smooth,'canny',TH);
% 
% imagesc(Edge_Map);
% 
% 
% %% To calculate the FF based on input image
% 
% 
% 
% 
% 
% 
% subplot(2,1,1);
% imagesc(Data_Image_Raw)
% axis equal
% colormap('gray')
% hold on
% plot([N,N],[size(Data_Image_Raw,2),1],'Color','r','LineWidth',2)
% 
% hold off
% %ylim([N-50 N+50])
% %xlim([400 550])
% 
% subplot(2,1,2);
% plot(Data_Image_Raw(:,N),'-*');
% 
% 
% %%
% for p=1:length(file_list)
%     file_name=file_list(p).name;
%     file_path=[folder_path file_name];
% 
%     Data_Save_Folder='D:\170328_Resolution test\Processed Data\';
%     mkdir(Data_Save_Folder);
%     Processed_Data_Path=[Data_Save_Folder Data_Name '_' file_name '.bin'];
%     cd(folder_path);
% 
%     Data_Image_Raw=mean(double(imread(file_path,'tiff')),3);
%     imagesc(Data_Image_Raw)
% end