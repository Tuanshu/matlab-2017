clear all

Phi_TiSa=20;
NA_TiSa=0.75;

Const=Phi_TiSa*asin(NA_TiSa);
B_TiSa=1/(Phi_TiSa*asin(NA_TiSa));

Phi_MMF=106;
NA_MMF=0.22;

Phi_MMF_Min=Const/asin(NA_MMF);
NA_MMF_Min=Const/Phi_MMF;
B_MMF=1/(Phi_MMF*asin(NA_MMF));


NA_Obj=0.5;
Phi_MMF_Obj_Min=1/B_MMF/asin(NA_Obj);

