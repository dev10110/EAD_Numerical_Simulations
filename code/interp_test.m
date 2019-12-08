clear;
clc;
close all;

load interp_test

x = linspace(-0.2,0.2,1000);
y = linspace(-0.1,0.2,1000);

tic
for i = 1:length(x)

    temp = uxinterp(x(i),y(i));

end
toc

tic
temp = uxinterp(x,y);
toc




