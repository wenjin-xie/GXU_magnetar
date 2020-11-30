
clc;
clear;
v=0.5;
w=1;
syms x;
r=real(v+w.*cos(x)*1i+w.*sin(x)*1i);
d=v+w.*cos(x)*1i+w.*sin(x)*1i;
d1=d/r
