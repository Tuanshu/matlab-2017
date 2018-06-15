clear all

f_c=30;

L_total_min=58.5+72.5-15;
L_total_max=58.5+72.5+15;

d1=31:0.01:100;

d2=1./((1/f_c)-(1./d1));



b2_4ac_min=(L_total_min^2-4*f_c*L_total_min)^0.5;   %??, ?????, ????4????!
b2_4ac_max=(L_total_max^2-4*f_c*L_total_max)^0.5;

d1_max_plus=0.5*(L_total_max+b2_4ac_max);
d1_max_minus=0.5*(L_total_max-b2_4ac_max);


%% to find the relation between d1 and L_total

L_total_Calculated=d2+d1;

plot(d1,L_total_Calculated);

[L_total_Calculated_min L_total_Calculated_minindex]=min(L_total_Calculated);

M=d2./d1;


plot(d1,L_total_Calculated);
hold on
plot(d1,M);
hold off

subplot(2,1,1)
plot(L_total_Calculated,M);
xlim([115 180]);
xlabel('Total Length Between Two Intermediate Plane (mm)');
ylabel('Magnification');

subplot(2,1,2)
plot(d2,M);
xlim([115 180]-f_c);    %Assume d1?f_c"??"
xlabel('d2 (mm)');
ylabel('Magnification');