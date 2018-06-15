clear all


fobj=9;
fproj=150;
d=150.0000;
Dtotal=160;
d_error=[0:1:100]; %mm
Line_Width=0.1; %mm

for p=1:length(d_error)
    pos_1=1/(1/fproj-1/(d+d_error(p)));

    d_2=Dtotal-pos_1;

    pos_2=1/(1/fobj-1/d_2);

    M1=pos_1/(d+d_error(p));
    M2=pos_2/d_2;

    M(p)=M1*M2
    GP1_Pos(p)=pos_2;
    Line_Width_Pixel(p)=abs(Line_Width/M(p)/7.4);
end


plot(d_error,1./M, 'LineWidth', 2);

xlabel('CCD displacement (mm)','fontsize',15);
ylabel('Magnification','fontsize',15);
set(gca,'fontsize',15)


plot(d_error,Line_Width_Pixel, 'LineWidth', 2);

xlabel('CCD displacement (mm)','fontsize',15);
ylabel('Magnification','fontsize',15);
set(gca,'fontsize',15)
