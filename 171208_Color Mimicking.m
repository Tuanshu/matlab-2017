clear all

Ori=imread('D:\000177.tiff')-45;
Modi=ColorMimicking(Ori,'D:\001368.tiff','tiff','tiff');


Modi_Image=Ori;

Modi_Image(:,:,1)=Modi_Image(:,:,1)*Modi(1);
Modi_Image(:,:,2)=Modi_Image(:,:,2)*Modi(2);
Modi_Image(:,:,3)=Modi_Image(:,:,3)*Modi(3);

imagesc(Modi_Image);
imwrite(Modi_Image,'D:\000177_2.tiff');
%%
Image=imread('O:\?????\???\171207_DiMirau Test\030-A_22fps_B130_R2.73_B1_No Sample_DiMirau_no BCG_10ms_Gain300_1\002010.tiff');

Mean=mean(Image(:))