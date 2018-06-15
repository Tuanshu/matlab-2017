clear all

R1_1=5;
R1_2=30;
R1_3=95;


R2_1=30;
R2_2=30;
R2_3=30;


A=1;
a=0.002;        %??tissue??
a_b=0.01;       %?tissue??
b=a_b-a;        %??tissue?????
c=0.007;        %GP2 stray light
d=0.775;        %selective coating?GP1??? (filling ratio)
Nsat=20000*16*4;    %因橫向上用4個frames平均

R1_Array=0:0.001:1;
R2_Array=0:0.001:1;

R1=repmat(R1_Array,[length(R2_Array) 1]);
R2=repmat(R2_Array',[1 length(R1_Array)]);

T_nonselect=((1-R1).^2).*(((1-R2).^2)*(a+b)+(R2.^2).*R1+R2*c);
T_select=d.^2*(((1-R2).^2)*(a+b)+(R2.^2).*R1+R2*c);

SNR=((A^2)*Nsat.*((1-R2).^2).*(R2.^2).*R1.*a)./((((1-R2).^2)*(a+b)+(R2.^2).*R1+R2*c).^2);

SNR_log=log10(SNR)*10;

imagesc(SNR_log,'xdata',R1_Array*100,'ydata',R2_Array*100)
set(gca,'YDir','normal')

xlabel('R1 (%)')
ylabel('R2 (%)')
caxis([0 max(SNR_log(:))])
axis equal
xlim([0 100])
ylim([0 100])
colorbar 

R1_index1=find(R1_Array*100>R1_1,1,'first');
R2_index1=find(R2_Array*100>R2_1,1,'first');
R1_index2=find(R1_Array*100>R1_2,1,'first');
R2_index2=find(R2_Array*100>R2_2,1,'first');
R1_index3=find(R1_Array*100>R1_3,1,'first');
R2_index3=find(R2_Array*100>R2_3,1,'first');
text(R1_1,R2_1,['．' sprintf('%2.2g dB',SNR_log(R2_index1,R1_index1))],'Color','white'); 
text(R1_2,R2_2,['．' sprintf('%2.2g dB',SNR_log(R2_index2,R1_index2))],'Color','white'); 
text(R1_3-10,R2_3,['．' sprintf('%2.2g dB',SNR_log(R2_index3,R1_index3))],'Color','white'); 


% 
% subplot(1,2,1)
% imagesc(T_nonselect,'xdata',R1_Array*100,'ydata',R2_Array*100)
% title('Non-selective coating Case')
% 
% set(gca,'YDir','normal')
% 
% xlabel('R1 (%)')
% ylabel('R2 (%)')
% caxis([0.01 0.1])
% axis equal
% xlim([0 100])
% ylim([0 100])
% colorbar 
% 
% text(R1_1,R2_1,['．' sprintf('T=%2.2g%%',T_nonselect(R2_index1,R1_index1)*100)],'Color','white'); 
% text(R1_2,R2_2,['．' sprintf('T=%2.2g%%',T_nonselect(R2_index2,R1_index2)*100)],'Color','white'); 
% 
% subplot(1,2,2)
% 
% imagesc(T_select,'xdata',R1_Array*100,'ydata',R2_Array*100)
% title('Selective coating Case')
% 
% set(gca,'YDir','normal')
% 
% xlabel('R1 (%)')
% ylabel('R2 (%)')
% caxis([0.01 0.1])
% axis equal
% xlim([0 100])
% ylim([0 100])
% colorbar 
% 
% text(R1_3-10,R2_3,['．' sprintf('T=%2.2g%%',T_select(R2_index3,R1_index3)*100)],'Color','white'); 