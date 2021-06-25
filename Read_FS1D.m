% Read ASCII output

clear all
close all
clc

% Open file
FileID = fopen('Gaussian_test-0000.dat');

% Read number of cells
IMAX = fscanf(FileID,'%s \n',1); 

% Read data
x   = fscanf(FileID,'%f \n',IMAX);
eta = fscanf(FileID,'%f \n',IMAX);
u   = fscanf(FileID,'%f \n',IMAX);

% Plot data
figure(1);
plot(x,eta,'o')
title('Free surface elevation')

figure(2);
plot(x,u,'o')
title('Velocity')
