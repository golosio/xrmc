clear all; close all
data=load('plot.dat');
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
for i=1:size(th,1)
	valm(i)=mean(val(i,:));
	sval(i)=std(val(i,:));
endfor
hth=th(10:17);
hval=zeros(8,1);
hvals=hval;
for i=1:8
	hval(9-i)=(valm(1+i)+valm(18-i))/2;
        hvals(9-i)=sqrt(sval(1+i)**2+sval(18-i)**2)/2;
endfor
semilogy(th1,valth,'r');
axis([0 90 5e-5 5e-3])
hold on
semilogyerr(hth,hval,hvals);
%semilogyerr(th(10:17),valm(10:17),sval(10:17));
figure(2)
semilogy(th1,valth,'r');
axis([-90 90 5e-5 5e-3])
axis([-100 100 5e-5 5e-3])
hold on
semilogyerr(th,valm,sval);
data2=[th valm sval];
%fp=fopen("valm.dat","w");
%fprintf(fp,"%e\t%e\t%e\n", data2');
%fclose(fp);
