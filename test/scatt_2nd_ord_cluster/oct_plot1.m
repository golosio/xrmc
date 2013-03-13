clear all; close all
data=load('plot1.dat');
th=data(:,1);
val=data(:,2:11);
valm=th;
sval=th;
for i=1:size(th,1)
	valm(i)=mean(val(i,:));
	sval(i)=std(val(i,:));
endfor
data1=load('theor_plot/theoretical.dat');
th1=data1(:,1);
valth=data1(:,2);
plot(th1,valth,'r');
axis([-100 100 0 1e-3])
hold on
errorbar(th(4:15),valm(4:15),sval(4:15));

