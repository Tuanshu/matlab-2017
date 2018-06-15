clear all
%%

Reference_Image_Path='D:\171116_Dermoscope Reference\IMG_0051_Test.jpg';
Sample_Image_Path='D:\171116_Post-processing of Scheme A Camera\Image0531.tiff';


Reference_Image=imread(Reference_Image_Path,'jpg');
Sample_Image=imread(Sample_Image_Path,'tiff');
Sample_Image=Sample_Image(:,:,1:3);

RGB_Ref=sum(sum(Reference_Image,1),2)./sum(sum(sum(Reference_Image,1),2),3);
RGB_Sam=sum(sum(Sample_Image,1),2)./sum(sum(sum(Sample_Image,1),2),3);

RGB_Modi=RGB_Ref./RGB_Sam;