clear all


fobj=9;
fproj=50;
d=50.00001;
D=60;
d_error=[-4:0.1:1]; %mm
CCD_Size=648*7.4;

for p=1:length(d_error)
    pos_1=1/(1/fproj-1/(d+d_error(p)));

    d_2=D-pos_1;

    pos_2=1/(1/fobj-1/d_2);

    M1=pos_1/(d+d_error(p));
    M2=pos_2/d_2;

    M=M1*M2
    GP1_Pos(p)=pos_2;
    Image_Size(p)=abs(CCD_Size*M);
end

plot(d_error,Image_Size-Image_Size(find(d_error>0,1,'first')), 'LineWidth', 2);

xlabel('CCD displacement (mm)','fontsize',15);
ylabel('Image Size Error (micron)','fontsize',15);
set(gca,'fontsize',15)

plot(d_error,Image_Size, 'LineWidth', 2);

xlabel('CCD displacement (mm)','fontsize',15);
ylabel('Image Size Error (micron)','fontsize',15);
set(gca,'fontsize',15)
%%

plot(d_error,(GP1_Pos-GP1_Pos(find(d_error>0,1,'first')))*1000, 'LineWidth', 2);

xlabel('CCD displacement (mm)','fontsize',15);
ylabel('GP1 displacement (micron)','fontsize',15);
set(gca,'fontsize',15)


plot(d_error,(GP1_Pos)*1000, 'LineWidth', 2);

xlabel('CCD displacement (mm)','fontsize',15);
ylabel('GP1 displacement (micron)','fontsize',15);
set(gca,'fontsize',15)