close all; clear all;
fp=fopen('output/output_1.dat','rb');
data=fread(fp, 800, 'float64');
fclose(fp);
I=data(401:800);
E=linspace(0,40,401);
E=E(2:401);
figure(1);
plot(E,I);

Sn_signal=zeros(41);
for i=0:40
   file_name=sprintf('output/output_%d.dat', i);
   file_name;
   fp=fopen(file_name,'rb');
   data=fread(fp, 800, 'float64');
   fclose(fp);
   Sn_signal(i+1)=sum(data(640:659));
end
figure(2);
plot(Sn_signal);

