clear all; close all;
data=load('output.dat');
I=data(2049:4096);
E=linspace(0,40.96,2048);
plot(E,I);
sum(I)
