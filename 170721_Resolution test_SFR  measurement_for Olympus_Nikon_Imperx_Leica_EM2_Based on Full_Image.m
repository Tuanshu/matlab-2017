clear all
%%
root_folder_path='C:\TuanShu\170721_Double Check SFR\tiff only\';

M=20;
Given_sub_folder_path='';

System='EM2PL_ExtendFOV';


ESF_Est_Type=2; % 1 and 2 for error function rought estimation, 3 for fermi
%root_folder_path='I:\161223_e2v EM1\CeYAG\Different Exp Time';
If_Best_Only=0;
Best_Index=65;

FermiFit = fittype( @(a, b, c, d, x) a./(exp(((x-b)/c))+1)+d);


Image_Downsampling_Ratio=1;
if strcmp(System,'Olympus')
    Camera_Pixel_Size=6.45; %micron
    Adaptor_M=0.63;
elseif strcmp(System,'NikonSystem-OlympusObj')
    Camera_Pixel_Size=3.45; %micron
    Adaptor_M=200/180;
elseif strcmp(System,'Imperx')
    Camera_Pixel_Size=7.4; %micron
    Adaptor_M=50/180;
elseif strcmp(System,'OlympusDP73')
    Camera_Pixel_Size=4.4; %micron
    Adaptor_M=0.5;
elseif strcmp(System,'Leica')
    Camera_Pixel_Size=6.5; %micron
    Adaptor_M=1;
elseif strcmp(System,'EM2Imperx')
    Camera_Pixel_Size=7.4; %micron
    Adaptor_M=150/180;
elseif strcmp(System,'EM2PL_ExtendFOV')
    Camera_Pixel_Size=10.6; %micron
    Adaptor_M=200/180;
end


Edge_TH=0.05;

Sampling_Resolution=Camera_Pixel_Size/M/Adaptor_M;
Edge_Window_Half_Size=round(50/Sampling_Resolution);


sub_processed_path='Processed Data';

Peak_Corr=0;

if strcmp(System,'EM2Imperx')
    SFR_Window_Half_Width=20;   %100
    Left_Bound=20;  %100
    Right_Bound=20; %100
elseif strcmp(System,'EM2PL_ExtendFOV')
    SFR_Window_Half_Width=20;   %100
    Left_Bound=20;  %100
    Right_Bound=20; %100
else
    SFR_Window_Half_Width=100;   %100
    Left_Bound=100;  %100
    Right_Bound=100; %100
end
Filter_Window_HalfWidth=10;
Further_Oversampling_Ratio=4;

%% If automatic sub folder
if strcmp(Given_sub_folder_path,'')

    Sub_Folder_List=dir(root_folder_path);
    
    for k = length(Sub_Folder_List):-1:1
        if ~Sub_Folder_List(k).isdir
            Sub_Folder_List(k) = [ ];
            continue
        end
        fname = Sub_Folder_List(k).name;
        if (fname(1) == '.') || strcmp(fname,'Processed Data') ||  strcmp(fname,'Other')
            Sub_Folder_List(k) = [ ];
        end
    end
    TotalFolderNumber=length(Sub_Folder_List);
else
    TotalFolderNumber=1;
end

for FolderNumber=1:TotalFolderNumber
    if strcmp(Given_sub_folder_path,'')
        sub_folder_path=Sub_Folder_List(FolderNumber).name;
        folder_path=[root_folder_path '\' sub_folder_path '\'];
    else
        sub_folder_path=Given_sub_folder_path;
        folder_path=[root_folder_path '\' sprintf('%s_%dx',sub_folder_path,M) '\'];

    end
    
    %folder_path=[root_folder_path '\' sprintf('%s_%dx',sub_folder_path,M) '\'];
    Data_Save_Folder=[root_folder_path '\' sub_processed_path '\'];
    if exist(Data_Save_Folder,'dir')==0
        mkdir(Data_Save_Folder);
    end
    
    if exist([Data_Save_Folder '\Resolution Chart\'],'dir')==0
        mkdir([Data_Save_Folder '\Resolution Chart\']);
    end

    Data_Name=sub_folder_path;

    %% Preliminary Data Reading
    file_list=dir(folder_path);

    for k = length(file_list):-1:1
        if file_list(k).isdir
            file_list(k) = [ ];
            continue
        end
        fname = file_list(k).name;
        if (fname(1) == '.') ||  strcmp(fname,'Thumbs.db')
            file_list(k) = [ ];
        end
    end

    if If_Best_Only == 1
        TotalFileNumber=1;
        Resolution_Best=99*ones(1,1);
        Resolution_Mean=99*ones(1,1);
        Resolution_STD=99*ones(1,1);
    else
        TotalFileNumber=floor(size(file_list,1)/Image_Downsampling_Ratio);
        %Resolution_Best=99*ones(size(file_list,1)/Image_Downsampling_Ratio,1);
        %因為可能multiple frame in one file, 所以算了
        %Resolution_Mean=99*ones(size(file_list,1)/Image_Downsampling_Ratio,1); 
        %Resolution_STD=99*ones(size(file_list,1)/Image_Downsampling_Ratio,1);
        Resolution_Best=[];
        %因為可能multiple frame in one file, 所以算了
        Resolution_Mean=[]; 
        Resolution_STD=[];
    end

    %%
    % Fermi Fit Type
    Acc_FrameNumber=1;
    for FileNumber=1:TotalFileNumber
        if If_Best_Only == 1
            file_name=file_list(Best_Index).name;
        else
            file_name=file_list(1+(FileNumber-1)*Image_Downsampling_Ratio).name;        
        end
        file_path=[folder_path file_name];
        TotalFrameNumber=length(imfinfo(file_path));
        for FrameNumber=1:TotalFrameNumber
            if strcmp(System,'Imperx')
                Data_Image_Raw=mean(double(imread(file_path,'bmp')),3);%.*FF;
                Data_Image_Raw=Data_Image_Raw(1:450,:);
            elseif strcmp(System,'EM2Imperx')
                Data_Image_Raw=double(imread(file_path,'tiff',FrameNumber))';%.*FF;
            elseif strcmp(System,'EM2PL_ExtendFOV')
                Data_Image_Raw=double(imread(file_path,'tiff',FrameNumber))';%.*FF;
                Data_Image_Raw=flip(Data_Image_Raw(101:900,486:555),1);
            else
                %Data_Image_Raw=mean(double(imread(file_path,'tiff')),3);%.*FF;
                Data_Image_Raw=double(imread(file_path,'tiff',FrameNumber));%.*FF;
            end
            % 2-bound edge detection

            X=1:size(Data_Image_Raw,1);

            L=mean(Data_Image_Raw(:,1:Left_Bound),2);
            R=mean(Data_Image_Raw(:,(size(Data_Image_Raw,2)-Right_Bound+1):size(Data_Image_Raw,2)),2);

            L_grad=[diff(L);0];
            R_grad=[diff(R);0];
            %plot(X,L,X,R,X,L_grad,X,R_grad);
            [L_Extreme_value L_Extreme]=findpeaks(abs(L_grad),'MINPEAKDISTANCE',round(80/Sampling_Resolution),'MINPEAKHEIGHT',max(abs(L_grad)/4));
            [R_Extreme_value R_Extreme]=findpeaks(abs(R_grad),'MINPEAKDISTANCE',round(80/Sampling_Resolution),'MINPEAKHEIGHT',max(abs(R_grad)/4));
            L_Extreme=L_Extreme+Peak_Corr;
            R_Extreme=R_Extreme+Peak_Corr;


            %% Find extremes with specific lateral spacing

            %% Take Gradient on binary gradient map
            X_grad_Filter=[-1*ones(Filter_Window_HalfWidth*2,Filter_Window_HalfWidth) 1*ones(Filter_Window_HalfWidth*2,Filter_Window_HalfWidth)];
            Y_grad_Filter=[-1*ones(Filter_Window_HalfWidth,Filter_Window_HalfWidth*2); 1*ones(Filter_Window_HalfWidth,Filter_Window_HalfWidth*2)];

            X_grad=filter2(X_grad_Filter,Data_Image_Raw,'valid');
            Y_grad=filter2(Y_grad_Filter,Data_Image_Raw,'valid');

            Abs_grad=(X_grad.^2+Y_grad.^2).^0.5;
            Abs_grad_sort=sort(Abs_grad(:),'descend');
            Grad_Top=mean(Abs_grad_sort(1:round(length(Abs_grad_sort)*Edge_TH)));

            Grad_Mask=Abs_grad>Grad_Top;

            Grad_Mask_Ske=bwmorph(Grad_Mask,'skel',Inf);
            Group_Edges=bwconncomp(Grad_Mask_Ske,8);
            Edge_Length_List=zeros(Group_Edges.NumObjects,1);
            for p=1:Group_Edges.NumObjects
                Edge_Length_List(p)=length(Group_Edges.PixelIdxList{p});
            end

            [maxlength maxlengthindex]=max(Edge_Length_List);
            Longest_Edge_1D=Group_Edges.PixelIdxList{maxlengthindex};
            [Longest_Edge_Y,Longest_Edge_X] = ind2sub([size(Grad_Mask_Ske,1) size(Grad_Mask_Ske,2)],Longest_Edge_1D);

            Matrix=[Longest_Edge_X ones(length(Longest_Edge_Y),1)];
            V=Matrix\Longest_Edge_Y;
            %Y_Fit=V(1)*Longest_Edge_X+V(2)
            Angle=atan(V(1))*180/pi;
            %Angle_Map=atan(X_grad./Y_grad)*180/pi;


            D_L_R=abs((1+Left_Bound)/2-(2*size(Data_Image_Raw,2)-Right_Bound+1)/2);
            LR_Extreme=L_Extreme+V(1)*D_L_R;

            %% 17/4/5 ??length(L_Extreme)==length(R_Extreme), ???????edge??(????), New Method: uese the closest
            %R_Extrme_First_Index=find(abs(LR_Extreme(1)-R_Extreme)<(600/M),1,'first');
            [value R_Extrme_First_Index]=min(abs(LR_Extreme(1)-R_Extreme));
            %L_Extrme_Last_Index=find(abs(LR_Extreme-R_Extreme(end))<(600/M),1,'last');
            [value L_Extrme_Last_Index]=min(abs(LR_Extreme-R_Extreme(end)));
            L_Extreme=L_Extreme(1:L_Extrme_Last_Index);
            R_Extreme=R_Extreme(R_Extrme_First_Index:end);

            Number_of_Edge=length(L_Extreme);

            if length(R_Extreme)>Number_of_Edge
                if abs(R_Extreme(1)-1)<abs(R_Extreme(end)-size(Data_Image_Raw,1))
                    R_Extreme(1)=[];
                else
                    R_Extreme(end)=[];
                end
            end


            if length(R_Extreme)~=Number_of_Edge
                warning('#of L&R Edge not Equal');
                continue;
            end

            %% Slope calculation (for SFR)
             Angle_Array=-1*atan((L_Extreme-R_Extreme)/D_L_R)*180/pi;
             Angle_Mean=mean(Angle_Array); %degree
             Angle_STD=std(Angle_Array);


            %% Cosine Correction to Sampling Frequency
            Sampling_Resolution_Corr=abs(Sampling_Resolution*cos(Angle_Mean/180*pi));

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


            %% Plot Edge on Image (real)
        % %     
        %     imagesc(Data_Image_Raw);
        %     colormap(gray);
        %     hold on
        %     for p=1:Number_of_Edge
        %         plot([(1+Left_Bound)/2,(2*size(Data_Image_Raw,2)-Right_Bound+1)/2],[L_Extreme(p),R_Extreme(p)],'Color','r','LineWidth',2)
        %     end
        %     plot([1:size(Y_Index_Array,2)],Y_Index_Array,'.');
        %     hold off
        %     saveas(gcf,[root_folder_path '\' sub_folder_path '.png'])

            %% Edge Window Defination (3D array for SFR, and DO NOT shift for different lateral position)

            Edge_Array=zeros(2*Edge_Window_Half_Size+1,SFR_Window_Half_Width*2+1,Number_of_Edge);

            for p=1:Number_of_Edge
                for q=1:(SFR_Window_Half_Width*2+1)
                        Edge_Array(:,q,p)=Data_Image_Raw((round(Y_Index_Array(p,size(Data_Image_Raw,2)/2))-Edge_Window_Half_Size):(round(Y_Index_Array(p,size(Data_Image_Raw,2)/2))+Edge_Window_Half_Size),size(Data_Image_Raw,2)/2-SFR_Window_Half_Width+q);
                end
            end

            %plot(Edge_Array(:,:,1)); %????????, ????pixel phase???(note??shift?, ?phase?>1 pixel)
            %% Generate pixel phase array
            Pixel_Phase_Per_Lateral_Pixel=1/tan((90-Angle_Mean)/180*pi);    %Note: unit: pixel)
            if abs(1/Pixel_Phase_Per_Lateral_Pixel)>(2*SFR_Window_Half_Width)
                disp('Not enough Pixels to complete a cycle.')
                continue;
            end
            Upsampling_Ratio=abs(floor(1/Pixel_Phase_Per_Lateral_Pixel))*Further_Oversampling_Ratio;
            Pixel_Phase_Array=repmat((0:Pixel_Phase_Per_Lateral_Pixel:Pixel_Phase_Per_Lateral_Pixel*((SFR_Window_Half_Width*2+1)-1))',[1 Number_of_Edge]);
            %% ??????, ?????interpolation
            Manual_Adj_Phase_Ratio=1;   %????1???edge??
            Based_Index_Array_Old=[1:(2*Edge_Window_Half_Size+1)];
            % Range Estimation (to avoid NaN)
            Max_Pixel_Phase_Array=Manual_Adj_Phase_Ratio*max(Pixel_Phase_Array(:));
            Min_Pixel_Phase_Array=Manual_Adj_Phase_Ratio*min(Pixel_Phase_Array(:));
            Min_Based_Index_Array_New=min(Based_Index_Array_Old-Min_Pixel_Phase_Array);
            Max_Based_Index_Array_New=max(Based_Index_Array_Old-Max_Pixel_Phase_Array);
            Based_Index_Array_New=Min_Based_Index_Array_New:1/Upsampling_Ratio:Max_Based_Index_Array_New;

            Edge_Array_Interp=zeros(length(Based_Index_Array_New),size(Edge_Array,2),size(Edge_Array,3));
            for p=1:Number_of_Edge
                for q=1:(SFR_Window_Half_Width*2+1)
                        Edge_Array_Interp(:,q,p)=interp1(Based_Index_Array_Old-Manual_Adj_Phase_Ratio*Pixel_Phase_Array(q,p),Edge_Array(:,q,p),Based_Index_Array_New); %????????, ????????
                end
            end
            %
            %plot(Edge_Array_Interp(:,:,1)); %?????extrap, ?NaN????
            %% Mean and Smooth

            Edge_Array_Interp_Mean=squeeze(mean(Edge_Array_Interp,2));

            % Edge Array Normalization (Normalization AFTER SFR)
            New_Sampling_Resolution=Sampling_Resolution_Corr/Upsampling_Ratio;
            Edge_Array_Normalized=Edge_Array_Interp_Mean;
            for p=1:size(Edge_Array_Interp_Mean,2)
                Edge_Array_Normalized(:,p)=(Edge_Array_Interp_Mean(:,p)-min(Edge_Array_Interp_Mean(:,p)))/(max(Edge_Array_Interp_Mean(:,p))-min(Edge_Array_Interp_Mean(:,p)));
            end
            %New_Position_Array=New_Sampling_Resolution*Based_Index_Array_New;
            New_Position_Array=Sampling_Resolution_Corr*Based_Index_Array_New;   %17/4/5 corrected

            %

            %plot(Edge_Array_Normalized(:,p))
            % ESF Resolution Calc
            Resolution_Array=99*ones(1,size(Edge_Array_Normalized,2));
            for p=1:size(Edge_Array_Normalized,2)
                if ESF_Est_Type == 1
                    if Edge_Array_Normalized(1,p)>Edge_Array_Normalized(end,p)
                        Resolution_Array(p)=(find(Edge_Array_Normalized(:,p)<0.1,1,'first')-find(Edge_Array_Normalized(:,p)<0.9,1,'first'))*New_Sampling_Resolution;
                    else
                        Resolution_Array(p)=(find(Edge_Array_Normalized(:,p)>0.9,1,'first')-find(Edge_Array_Normalized(:,p)>0.1,1,'first'))*New_Sampling_Resolution;
                    end
                elseif ESF_Est_Type == 2
                    if Edge_Array_Normalized(1,p)>Edge_Array_Normalized(end,p)
                        Resolution_Array(p)=(find(Edge_Array_Normalized(:,p)<0.282,1,'first')-find(Edge_Array_Normalized(:,p)<0.7259,1,'first'))*New_Sampling_Resolution*2;
                    else
                        Resolution_Array(p)=(find(Edge_Array_Normalized(:,p)>0.7259,1,'first')-find(Edge_Array_Normalized(:,p)>0.282,1,'first'))*New_Sampling_Resolution*2;
                    end            
                elseif ESF_Est_Type == 3
                    if Edge_Array_Normalized(1,p)>Edge_Array_Normalized(end,p)
                        FermiResult=fit(New_Position_Array',Edge_Array_Normalized(:,p),FermiFit,'StartPoint', [1, 50, 100, 0]);
                        Resolution_Array(p)=abs(FermiResult.c)*3.5254;
                    else
                        FermiResult=fit(New_Position_Array',Edge_Array_Normalized(size(Edge_Array_Normalized,1):(-1):1,p),FermiFit, 'StartPoint', [1, 50, 100, 0]);
                        Resolution_Array(p)=abs(FermiResult.c)*3.5254;
                    end               
                end

            end
            min(Resolution_Array)

            %% Fermi Fitting
        %     FermiCurve=FermiResult.a./(exp((New_Position_Array'-FermiResult.b)/FermiResult.c)+1)+FermiResult.d;
        %     plot(New_Position_Array',FermiCurve,New_Position_Array',Edge_Array_Normalized(size(Edge_Array_Normalized,1):(-1):1,p));
        %     plot(New_Position_Array',FermiCurve,New_Position_Array',Edge_Array_Normalized(:,p));
        % 
        %     ylim([min(Edge_Array_Normalized(size(Edge_Array_Normalized,1):(-1):1,p)) max(Edge_Array_Normalized(size(Edge_Array_Normalized,1):(-1):1,p))])
        %     xlim([min(New_Position_Array) max(New_Position_Array) ]);

            %%
            % for p=1:size(Edge_Array_Normalized,2)
            %     plot(Edge_Array_Normalized(:,p));
            %     hold on
            % end
            plot(New_Position_Array,Edge_Array_Normalized);
            xlim([New_Position_Array(1) New_Position_Array(end)]);
            saveas(gcf,[root_folder_path '\' sub_processed_path '\' sub_folder_path sprintf('_Fig%d',Acc_FrameNumber) '_Butterfly.png'])
            Resolution_Best(Acc_FrameNumber)=min(Resolution_Array);
            Resolution_Mean(Acc_FrameNumber)=mean(Resolution_Array);
            Resolution_STD(Acc_FrameNumber)=std(Resolution_Array);
            disp(FrameNumber);
            Acc_FrameNumber=Acc_FrameNumber+1;
        end
    end

    if If_Best_Only == 1
        if ESF_Est_Type == 3
                FermiCurve=FermiResult.a./(exp((New_Position_Array'-FermiResult.b)/FermiResult.c)+1)+FermiResult.d;
                if Edge_Array_Normalized(1,p)>Edge_Array_Normalized(end,p)
                    plot(New_Position_Array',FermiCurve,New_Position_Array',Edge_Array_Normalized(:,p));
                else
                    plot(New_Position_Array',FermiCurve,New_Position_Array',Edge_Array_Normalized(size(Edge_Array_Normalized,1):(-1):1,p));
                end
                ylim([min(Edge_Array_Normalized(size(Edge_Array_Normalized,1):(-1):1,p)) max(Edge_Array_Normalized(size(Edge_Array_Normalized,1):(-1):1,p))])
                xlim([min(New_Position_Array) max(New_Position_Array) ]);
        end
        Resolution_Best
    else
        plot(Resolution_Best)
        xlabel('FileNumber');
        ylabel('Resolution (micron)');
        if ESF_Est_Type == 3
            saveas(gcf,[root_folder_path '\' sub_processed_path '\Resolution Chart\Resolution Chart_' sub_folder_path sprintf('_Min_%g_Fermi.png',min(Resolution_Best))])

        else
            saveas(gcf,[root_folder_path '\' sub_processed_path '\Resolution Chart\Resolution Chart_' sub_folder_path sprintf('_Min_%g.png',min(Resolution_Best))])
        end
    end  
end

    % %% Find The Best
% Smooth_Window=10;
% Resolution_Mean_Smooth=smooth(Resolution_Mean,Smooth_Window);
% 
% plot(Resolution_Mean_Smooth)
% xlabel('FileNumber');
% ylabel('Resolution (micron)');
% 
% [Value_of_the_best_resolution Index_of_the_best_resolution]=min(Resolution_Mean_Smooth);
% 
% %% Serious calc at the best
% if If_Best_Only ==1
%     
% end
% 
% 
% %%
%  
%  NN=14;
%  
%  plot(Edge_Array_Normalized(:,NN));
