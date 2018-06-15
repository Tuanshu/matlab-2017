clear all
%%
root_folder_path='C:\TuanShu\170613_PhotonFocus PTC Data Used\Used\';
Data_Name='170613_PhotonFocus_PTC';
Current_Array=[0.1:0.1:0.7];


Size=25;

H_Offset=300;
H_Range=30;
X_ROI=[H_Offset H_Offset+H_Range];


V_Offset=1;
V_Range=7;
Y_ROI=[V_Offset V_Offset+V_Range];

Noise_Removal_Method=4; %1: do nothing, 2: normalization, 3: filtering, 4: 4-point

% Feature for Filtering mode
Filtering_Ratio=0.02;       %意思是, 若Array長度為A, 則將外側(低頻)之A*Filtering_Ratio個pixels設為0
% Feature for 4-point mode
Averaging_Factor=1;
N=4;          %for N-point
Unbias_Estimator=0.9213;

Mean_Array_Map=zeros(X_ROI(2)-X_ROI(1)+1,Y_ROI(2)-Y_ROI(1)+1,length(Current_Array));
STD_Array_Map=zeros(X_ROI(2)-X_ROI(1)+1,Y_ROI(2)-Y_ROI(1)+1,length(Current_Array));

Mean_Array_for_norm=zeros(1,length(Current_Array));

Row=1024;
Colomn=1024;

for p=1:length(Current_Array)
    last_folder_name=sprintf('it%g_1',Current_Array(p));
    folder_path=[root_folder_path '\' last_folder_name];

    Data_Save_Folder='C:\TuanShu\Processed Data\';

    cd(folder_path);
    file_list=dir(folder_path);
    Original_Lnegth_file_list=length(file_list);
    for q=1:Original_Lnegth_file_list
        if file_list(Original_Lnegth_file_list+1-q).isdir ~= 0
            file_list(Original_Lnegth_file_list+1-q)=[];
        elseif strcmp(file_list(Original_Lnegth_file_list+1-q).name,'.') == 1
            file_list(Original_Lnegth_file_list+1-q)=[];
        elseif strcmp(file_list(Original_Lnegth_file_list+1-q).name,'..') == 1
            file_list(Original_Lnegth_file_list+1-q)=[];
        elseif strcmp(file_list(Original_Lnegth_file_list+1-q).name,'Processed Data') == 1
            file_list(Original_Lnegth_file_list+1-q)=[];
        end
    end
    Image_Stack=zeros(Row,Colomn,length(file_list));
    for r=1:length(file_list)
        
        file_path=[folder_path '\' file_list(r).name];
        fin=fopen(file_path);
        Image_Stack(:,:,r)=fread(fin,[Row,Colomn],'uint16'); %*Frame   不知為何, 看起來就像是要除16
        %%
        %fclose(fin);
        fclose('all');
    end
    

    Number_of_Filter_Pixel_one_side=size(Image_Stack,3)*Filtering_Ratio/2;   %總共有左右兩個side
    
    Mean_Array_for_norm=mean(mean(Image_Stack,2),1); % 固定用全張影像做, 有好有壞, 優點: 不會因為sample點數太少而不準, 缺點: 若影像有部分飽和則不準
    All_Mean=mean(Mean_Array_for_norm);
    Mean_Stack_for_norm=repmat(Mean_Array_for_norm,[size(Image_Stack,1) size(Image_Stack,2) 1]);
    Image_Stack_norm=Image_Stack./Mean_Stack_for_norm*All_Mean;
    
    Image_Stack_FFT=fft((Image_Stack-Mean_Stack_for_norm),[],3);
    
    Image_Stack_FFT_Filter=Image_Stack_FFT;
    Image_Stack_FFT_Filter(:,:,1:Number_of_Filter_Pixel_one_side)=0;
    Image_Stack_FFT_Filter(:,:,(size(Image_Stack_FFT_Filter,1)-Number_of_Filter_Pixel_one_side+1):size(Image_Stack_FFT_Filter,1))=0;
    Image_Stack_Filter=real(ifft(Image_Stack_FFT_Filter,[],3))+Mean_Stack_for_norm;

    Averaged_Length=floor(size(Image_Stack,3)/Averaging_Factor);
    Temp=0;
    for q=1:Averaging_Factor
       Temp=Temp+Image_Stack(:,:,(Averaging_Factor-(q-1)):Averaging_Factor:(Averaging_Factor*Averaged_Length)-(q-1));
    end
    Image_Stack_Ave=Temp/Averaging_Factor;
    Image_Stack_Npoint=zeros(size(Image_Stack_Ave,1),size(Image_Stack_Ave,2),floor(size(Image_Stack_Ave,3)/N));
    for q=1:size(Image_Stack_Npoint,3)
        Temp=Image_Stack_Ave(:,:,((q-1)*N+1):(q*N));
        Image_Stack_Npoint(:,:,q)=((N*sum(Temp.^2,3)-sum(Temp,3).^2).^0.5)*(2^0.5)/N;
    end
    
    
    if Noise_Removal_Method == 1
        Mean_Array_Map(:,:,p)=mean(Image_Stack(X_ROI(1):X_ROI(2),Y_ROI(1):Y_ROI(2),:),3);
        STD_Array_Map(:,:,p)=std(Image_Stack(X_ROI(1):X_ROI(2),Y_ROI(1):Y_ROI(2),:),0,3);
    elseif Noise_Removal_Method == 2
        Mean_Array_Map(:,:,p)=mean(Image_Stack_norm(X_ROI(1):X_ROI(2),Y_ROI(1):Y_ROI(2),:),3);
        STD_Array_Map(:,:,p)=std(Image_Stack_norm(X_ROI(1):X_ROI(2),Y_ROI(1):Y_ROI(2),:),0,3);    
    elseif Noise_Removal_Method == 3
        Mean_Array_Map(:,:,p)=mean(Image_Stack_Filter(X_ROI(1):X_ROI(2),Y_ROI(1):Y_ROI(2),:),3);
        STD_Array_Map(:,:,p)=std(Image_Stack_Filter(X_ROI(1):X_ROI(2),Y_ROI(1):Y_ROI(2),:),0,3);        
    elseif Noise_Removal_Method == 4
        Mean_Array_Map(:,:,p)=mean(Image_Stack(X_ROI(1):X_ROI(2),Y_ROI(1):Y_ROI(2),:),3);
        STD_Array_Map(:,:,p)=mean(Image_Stack_Npoint(X_ROI(1):X_ROI(2),Y_ROI(1):Y_ROI(2),:),3)*Unbias_Estimator;        
    end
    disp(p);
end
All_Meam=mean(Mean_Array_for_norm);

% %%
% NN=6;
% 
% plot(real(Data_Image_FFT(:,NN)));
% 

%


Mean_Array_Map_1D=Mean_Array_Map(:);
STD_Array_Map_1D=STD_Array_Map(:);
VAR_Array_Map_1D=STD_Array_Map_1D.^2;
scatter(Mean_Array_Map_1D,VAR_Array_Map_1D,Size,'filled');
xlabel('Signal (DN)','fontsize',15);
ylabel('Variance (DN^2)','fontsize',15);
set(gca,'fontsize',15)
ylim([0 200]);
xlim([0 4096]);
saveas(gcf,[Data_Save_Folder Data_Name sprintf('_Method_%g',Noise_Removal_Method) '.png']);
