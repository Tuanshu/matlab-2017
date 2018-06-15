clear all;

Number_of_Frame=400;
Header_Size=1024*1024;
Row_per_FOV=648;
Colume_per_FOV=488;

file_path='F:\20161230\1430548_80um\FOV Data\20161229193511\017_014.bin';
fin=fopen(file_path);
fseek(fin, Header_Size, 'bof');
fov=zeros(Row_per_FOV,Colume_per_FOV,Number_of_Frame);
for p=1:Number_of_Frame
    fov(:,:,p)=(fread(fin,[Row_per_FOV Colume_per_FOV],'single'));
end
Max=max(fov(:))
Min=min(fov(:))
Mean=mean(fov(:))
Std=std(fov(:))
a=hist(fov(:),Max)

a_total=sum(a)
a_500=sum(a(1:500))
b = cumsum(a)/a_total
plot(b,'LineWidth',3)
xlabel('DN','FontSize',12)
ylabel('H','FontSize',12)
set(gca,'FontSize',12);

ylabel('C','FontSize',12)
