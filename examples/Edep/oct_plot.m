clear all; close all;
data=load('profile_ph.txt');
data1=load('profile_ph_th.txt');
y=data(:,2);
y1=data1(:,2);
plot(y)
hold on
figure(1);
plot(y1,'r')

data2=load('profile_C.txt');
data3=load('profile_C_th.txt');
y2=data2(:,2);
y3=data3(:,2);
figure(2);
plot(y2)
hold on
plot(y3,'r')

data4=load('profile_tot.txt');
y4=data4(:,2);
figure(3);
plot(y4)
hold on
plot(y1+y3,'r')
