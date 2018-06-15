clear all
%%
If_Selective=1;
R1=1;
R2=0:0.01:1;
R_sam_Total=0.01;
Sam_Stary_Ratio=194;
R_sam=R_sam_Total/(Sam_Stary_Ratio+1);
R_sam_stray=R_sam_Total/(Sam_Stary_Ratio+1)*Sam_Stary_Ratio;
R_Ref_stray=0;
R_Constant_stray=0;

Ratio_Ref=R2.^2;

Ratio_Sam=(1-R2).^2;

P_Ref=R1*Ratio_Ref;
P_Sam=R_sam*Ratio_Sam;


P_Inter=2*(P_Ref.*P_Sam).^0.5;

P_Sam_Stray=R_sam_stray.*Ratio_Sam;
P_Ref_Stray=R_Ref_stray.*Ratio_Ref;

IE=P_Inter./(P_Ref+P_Sam+P_Sam_Stray+P_Ref_Stray+R_Constant_stray);

P_Sam_Total=P_Sam+P_Sam_Stray;
P_Ref_Total=P_Ref+P_Ref_Stray;

plot(R2,IE,R2,P_Sam_Total,R2,P_Ref_Total);
plot(R2,IE);