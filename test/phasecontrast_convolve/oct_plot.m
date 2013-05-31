close all; clear all;
data1=load('profilo.dat');
x1=(data1(:,1)-400+0.5)*5e-5;
I1=data1(:,2);
plot(x1, I1)
hold on
data2=load('arfelli.dat');
x2=data2(:,1);
I2=data2(:,2)*4.23e-10;
plot(x2, I2,'r')

