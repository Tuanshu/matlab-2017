clear all

    H_Offset_View=190;
    H_Range_View=650;
    Horizontal_ROI_View=[H_Offset_View H_Offset_View+H_Range_View];

    %
    H_Offset=360;
    H_Range=460;
    Horizontal_ROI=[H_Offset H_Offset+H_Range];

    Data_Image_Raw=double(imread('I:\170118_TiO2 Test\Water_Again_14\000009.tiff','tiff'));

    Data_Image_Raw_Cut=Data_Image_Raw(:,Horizontal_ROI_View(1):Horizontal_ROI_View(2));
    
    Data_Image_Raw_Norm=Data_Image_Raw_Cut./repmat(mean(Data_Image_Raw(:,Horizontal_ROI(1):Horizontal_ROI(2)),2),[1 size(Data_Image_Raw_Cut,2)]);
    
    Data_Image_Raw_BND=repmat(mean(Data_Image_Raw_Norm,1),[10000 1]);
    
    Data_Image_Raw_Sub=Data_Image_Raw_Norm-Data_Image_Raw_BND;
    
    imagesc(Data_Image_Raw_Sub);
    colormap(gray);
%% STFT
    N_of_Line=577;
    Window=1000;
    
    Line_Profile=Data_Image_Raw_Sub(:,N_of_Line);
    
    X=spectrogram(Line_Profile,Window);