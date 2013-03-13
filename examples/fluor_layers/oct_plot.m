clear all; close all;
data=load('convolved_spectrum.dat');
E=data(:,1);
I=data(:,2);
semilogy(E,I);
hold on
data1=load('experimental.dat');
E1=data1(:,1);
I1=data1(:,2);
semilogy(E1,I1*2e-5,'r');
axis([2 30 1e-4 1e2]);
