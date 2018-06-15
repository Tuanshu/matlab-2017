clear all;

Sample=3;

if Sample==1
    folder_path='I:\20161025_BCCfresh_2\A1\BC23_0_3.2\FOV Data\20161025123800\';
    save_folder_path='I:\20161025_BCCfresh_2\A1\BC23_0_3.2\';
    start_number_of_row=6;
    start_number_of_colume=8;
    end_number_of_row=17;
    end_number_of_colume=23;
elseif Sample==2
    folder_path='I:\20161025_BCCfresh_2\A2\CD23_-10.2_-7_good\FOV Data\20161025121419\';
    save_folder_path='I:\20161025_BCCfresh_2\A2\CD23_-10.2_-7_good\';
    start_number_of_row=12;
    start_number_of_colume=8;
    end_number_of_row=23;
    end_number_of_colume=23;
elseif Sample==3
    folder_path='I:\20161025_BCCfresh_2\A3\BC23_-3.2_0\FOV Data\20161025122416\';
    save_folder_path='I:\20161025_BCCfresh_2\A3\BC23_-3.2_0\';
    start_number_of_row=6;
    start_number_of_colume=8;
    end_number_of_row=17;
    end_number_of_colume=23;
end
Header_Size=1;
num_of_row=end_number_of_row-start_number_of_row+1;
num_of_colume=end_number_of_colume-start_number_of_colume+1;
Row_per_FOV=648;
Colume_per_FOV=488;
row_pixel=num_of_row*Row_per_FOV;
colume_pixel=num_of_row*Colume_per_FOV;
frame_width=Row_per_FOV;
frame_height=Colume_per_FOV;

%%  read exposure time
exposure=zeros(num_of_row,num_of_colume);
for i=start_number_of_row:end_number_of_row
    for j=start_number_of_colume:end_number_of_colume
        file_path=[folder_path sprintf('%03d_%03d.bin',i,j)];
        fin=fopen(file_path);
        Skip_length=160;
        fseek(fin, Skip_length, 'bof');
        Header=fread(fin,Header_Size,'ulong');
        exposure((i-start_number_of_row+1),(j-start_number_of_colume+1))=Header;
        fclose(fin);
        disp(i*num_of_row+j);
    end
end
figure
imagesc(exposure)

exposuremap=zeros(row_pixel,colume_pixel);
for i=1:num_of_row
    for j=1:num_of_colume
        ff=exposure(i,j);
        exposuremap((648*(i-1)+1):648*i,(488*(j-1)+1):488*j)=ff;
        disp(i*num_of_row+j);
    end
end


%%  read fov data
skip_frame=0;
num_of_frame=16;

FOV_Record=zeros(frame_width,frame_height,num_of_row*num_of_colume);

row_Array=start_number_of_row:end_number_of_row;
column_Array=start_number_of_colume:end_number_of_colume;

for i=1:num_of_row
    for j=1:num_of_colume
        Header_Size=1048576+Row_per_FOV*Colume_per_FOV*4*skip_frame;
        fov1=zeros(Colume_per_FOV,Row_per_FOV);
        for k=1:num_of_frame
            file_path=[folder_path sprintf('%03d_%03d.bin',row_Array(i),column_Array(j))];
            fin=fopen(file_path);
            fseek(fin, Header_Size, 'bof');
            %fov=fread(fin,[Row_per_FOV,Colume_per_FOV*num_of_frame],'float32','b');
            fov=(fread(fin,[Row_per_FOV Colume_per_FOV],'single'))';
            fclose(fin);
            Header_Size=Header_Size+Row_per_FOV*Colume_per_FOV*4;
            fov1=fov1+1/num_of_frame*fov;
        end
        fov2=imrotate(fov1,90);
        fov3=fliplr(fov2);
        %% Record Frame
        FOV_Record(:,:,((i-1)*num_of_colume)+j)=fov3;
        disp((i-1)*num_of_colume+j);

    end
end


%%
QQ=36;
imagesc(FOV_Record(:,:,QQ));
%%
TH=8;
TH_map=ones(Row_per_FOV,Colume_per_FOV)*TH;
DF_map=max(min(FOV_Record,[],3),TH_map);
    imagesc(DF_map);

        colormap('gray');
        %caxis([0 1]);
        axis equal
        xlim([0 size(DF_map,2)]);
        ylim([0 size(DF_map,1)]);
%%  DF 

    Times=4;
    Search_Window_Half_Size=10;
    Search_Window_Weight=fspecial('gaussian', Search_Window_Half_Size*2+1,Search_Window_Half_Size/2);    %gaussian mask
    Search_Window_Weight=Search_Window_Weight/max(Search_Window_Weight(:)); %因為是要取max而非平均, 所以filter的最大值應該是1
    Dark_frame_Temp=DF_map;
    Dark_frame_Out=DF_map;
    for w=1:Times
        for p=1:size(Dark_frame_Temp,1)
            for q=1:size(Dark_frame_Temp,2)
                x_min=max((p-Search_Window_Half_Size+1),1);
                x_max=min((p+Search_Window_Half_Size),size(Dark_frame_Temp,1));
                y_min=max((q-Search_Window_Half_Size+1),1);
                y_max=min((q+Search_Window_Half_Size),size(Dark_frame_Temp,2));
                Window_Temp=Search_Window_Weight((x_min-p+Search_Window_Half_Size+1):(x_max-(p+Search_Window_Half_Size)+Search_Window_Half_Size*2+1),(y_min-q+Search_Window_Half_Size+1):(y_max-(q+Search_Window_Half_Size)+Search_Window_Half_Size*2+1)).*Dark_frame_Temp(x_min:x_max,y_min:y_max);
                Dark_frame_Out(p,q)=max(Window_Temp(:));
            end
        end
        Dark_frame_Temp=Dark_frame_Out;
        disp(w);
    end
    
 Dark_frame=Dark_frame_Temp;
    
 DF_Reshape=reshape(Dark_frame,[size(Dark_frame,1)*size(Dark_frame,2) 1]);
 fid = fopen([save_folder_path 'DF.bin'], 'a+');
 fwrite(fid, DF_Reshape, 'single');
 fclose(fid);
 dlmwrite([save_folder_path,'Dark_frame.txt'],Dark_frame,'delimiter','\t','newline','pc','precision', '%.6f');
 
    subplot(1,2,1)

    imagesc(DF_map);
        colormap('gray');
        %caxis([0 1]);
        axis equal
        xlim([0 size(DF_map,2)]);
        ylim([0 size(DF_map,1)]);

    subplot(1,2,2)

    imagesc(Dark_frame);
        colormap('gray');
        %caxis([0 1]);
        axis equal
        xlim([0 size(Dark_frame_Out,2)]);
        ylim([0 size(Dark_frame_Out,1)]);






%% FF
Total_Image=zeros(Row_per_FOV,Colume_per_FOV);
Total_Number_of_Frame=num_of_row*num_of_colume;

Sub_Order=2;
Factor=1.2;

DF_Volume=repmat(Dark_frame,[1 1 Total_Number_of_Frame]);
if Sub_Order == 1
    FOV_DF_Sub=FOV_Record-DF_Volume*Factor;
elseif Sub_Order == 2
    FOV_DF_Sub_2=FOV_Record.^2-(DF_Volume*Factor).^2;
    FOV_DF_Sub_2(FOV_DF_Sub_2<0)=0;
    FOV_DF_Sub=FOV_DF_Sub_2.^0.5;
    FOV_DF_Sub(isnan(FOV_DF_Sub))=0;
end

%%%%%%
Averaged_Image=mean(FOV_DF_Sub,3);
FF_map=1./(Averaged_Image/mean(Averaged_Image(:)));
%%%%%%%

imagesc(Averaged_Image);
axis equal
fclose('all');
xlim([0 Colume_per_FOV]);
ylim([0 Row_per_FOV]);
colormap(gray);
axis off

imagesc(FF_map);
axis equal
fclose('all');
xlim([0 Colume_per_FOV]);
ylim([0 Row_per_FOV]);
colormap(gray);
axis off


FF_Reshape=reshape(FF_map,[size(FF_map,1)*size(FF_map,2) 1]);

 fid = fopen([save_folder_path 'FF.bin'], 'a+');
 fwrite(fid, FF_Reshape, 'single');
 fclose(fid);
 dlmwrite([save_folder_path,'FF.txt'],FF_map,'delimiter','\t','newline','pc','precision','%.6f');


 
 
 

 %% Overlapping Related

X_overlapping=25;
Y_overlapping=21;
frame_width_eff=frame_width-X_overlapping;
frame_height_eff=frame_height-Y_overlapping;

correction_A=ones(frame_width,frame_height);

for tt=1:X_overlapping
    correction_A(tt,:)=correction_A(tt,:)*((tt-1)/((X_overlapping-1)));
    correction_A(frame_width-tt+1,:)=correction_A(frame_width-tt+1,:)*((tt-1)/((X_overlapping-2))); 

end
for tt=1:Y_overlapping
    correction_A(:,tt)=correction_A(:,tt)*((tt-1)/(Y_overlapping-1));
    correction_A(:,frame_height-tt+1)=correction_A(:,frame_height-tt+1)*((tt-1)/((Y_overlapping-2))); 
end
 %% M
%M=zeros(frame_width_eff*num_of_row+X_overlapping,frame_height_eff*num_of_colume+Y_overlapping);
Image_Correlcted=zeros(frame_width_eff*num_of_row+X_overlapping,frame_height_eff*num_of_colume+Y_overlapping);


row_Array=start_number_of_row:end_number_of_row;
column_Array=start_number_of_colume:end_number_of_colume;


for i=1:num_of_row
    for j=1:num_of_colume
        if Sub_Order == 1
            Current=FOV_Record(:,:,((i-1)*num_of_colume)+j)-Dark_frame*Factor;
            Image_Correlcted(((i-1)*frame_width_eff+1):((i-1)*frame_width_eff+frame_width),((j-1)*frame_height_eff+1):((j-1)*frame_height_eff+frame_height))=Image_Correlcted(((i-1)*frame_width_eff+1):((i-1)*frame_width_eff+frame_width),((j-1)*frame_height_eff+1):((j-1)*frame_height_eff+frame_height))+Current.*FF_map.*correction_A;
        elseif Sub_Order == 2
            Current_2=FOV_Record(:,:,((i-1)*num_of_colume)+j).^2-(Dark_frame*Factor).^2;
            Current_2(Current_2<0)=0;
            Current=Current_2.^0.5;
            Image_Correlcted(((i-1)*frame_width_eff+1):((i-1)*frame_width_eff+frame_width),((j-1)*frame_height_eff+1):((j-1)*frame_height_eff+frame_height))=Image_Correlcted(((i-1)*frame_width_eff+1):((i-1)*frame_width_eff+frame_width),((j-1)*frame_height_eff+1):((j-1)*frame_height_eff+frame_height))+Current.*FF_map.*correction_A;
        end
    end
end
imagesc(Image_Correlcted)
caxis([5 25]);
colormap gray

%%
Min=5;
Max=25;    %因目前沒有平均, 所以訊號越強, 範圍就要越大

Image_Correlcted_Norm=(Image_Correlcted-Min)/((Max-Min));
Image_Correlcted_Norm(Image_Correlcted_Norm>1)=1;

Image_Correlcted_Norm(Image_Correlcted_Norm<0)=0;


imagesc(Image_Correlcted_Norm)
colormap gray

% multiply exposure
mkdir([sprintf('%s\\',save_folder_path),'processing']);
imwrite(Image_Correlcted_Norm,[sprintf('%s\\',save_folder_path),sprintf('processing\\processed_image_%d_%d_overlapping.png',X_overlapping,Y_overlapping)]);