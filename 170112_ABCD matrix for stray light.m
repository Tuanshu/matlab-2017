clear all

fobj=9;
d_GP2=7.35;
D=50;
fproj=150;
d_camera=fproj;

%M_x=[1 x; 0 1];                     %% MMF image to principal plane of WIO
M_obj=[1 0; -1/fobj 1];             %% WIO
M_GP2=[1 d_GP2; 0 1];               %% distance between WIO and GP2
M_GP1=[1 fobj; 0 1];               %% distance between WIO and GP2
M_D=[1 D; 0 1];                     %% distance between WIO and projection lens
M_proj=[1 0; -1/fproj 1];           %% projection lens
M_camera=[1 d_camera; 0 1];                %% distance between projection lens and camera

%M_stray_light=M_camera*M_proj*M_D*M_obj*M_GP2*M_GP2*M_obj*M_x;
M_stray_light_withoutX=M_camera*M_proj*M_D*M_obj*M_GP2*M_GP2*M_obj
M_stray_light_afterObj=M_obj*M_GP2*M_GP2*M_obj

% Assume x is known

x=-15.54546;            %-15.54546
M_x=[1 x; 0 1];                     %% MMF image to principal plane of WIO
Secondary_Image_Lateral_Offset=1;
Secondary_Image_Divergence=0.02778; %% 如果我要half width of FOI = 0.25
Input_Ray=[Secondary_Image_Lateral_Offset;Secondary_Image_Divergence];
M_stray_light=M_camera*M_proj*M_D*M_obj*M_GP2*M_GP2*M_obj*M_x;
Output_Ray=M_stray_light*Input_Ray
%%
% Calculate the corresponding FOI position (on both object and image plane) accroding to Secondary_Image_Lateral_Offset
%由於Secondary_Image_Divergence和FOI之寬度有關, 故或許可從這裡反推Secondary_Image_Divergence
Object_Ray=M_GP1*M_obj*M_x*Input_Ray        %1st element為half width of FOI
%Object_Ray_beforeProgress=M_obj*M_x*Input_Ray
%%
Marginal_Ray_afterX=

fobj=9;
fproj=50;
d=50;
D=0;
d_error=[50];

    for p=1:length(d_error)
    Image_Size=648*7.4; %micron

    Light_Image=[0;1];

    M_d=[1 d; 0 1];
    M_obj=[1 0; -1/fobj 1];
    M_proj=[1 0; -1/fproj 1];
    M_D=[1 D; 0 1];
    M_d_plus=[1 d+d_error(p)/2; 0 1];
    M_d_minus=[1 d-d_error(p)/2; 0 1];

    Light_Obj=M_obj*M_D*M_proj*M_d*Light_Image;
    Light_Obj_error_plus=M_obj*M_D*M_proj*M_d_plus*Light_Image;
    Light_Obj_error_minus=M_obj*M_D*M_proj*M_d_minus*Light_Image;

    Obj_Pos=Light_Obj(1)/Light_Obj(2)
    Obj_Pos_error_plus=Light_Obj_error_plus(1)/Light_Obj_error_plus(2)
    Obj_Pos_error_minus=Light_Obj_error_minus(1)/Light_Obj_error_minus(2)


    M=abs(Light_Obj(2)/Light_Image(2));

    M_error_plus=abs(Light_Obj_error_plus(2)/Light_Image(2));
    M_error_minus=abs(Light_Obj_error_minus(2)/Light_Image(2));

    Ratio=abs(M_error_plus-M_error_minus)/M;

    Size_Error(p)=abs(Image_Size/M_error_plus-Image_Size/M_error_minus);
    end

    plot(d_error,Size_Error);
    xlabel('CCD tolerance (mm)');
    ylabel('Object Space Size Error (micron)');