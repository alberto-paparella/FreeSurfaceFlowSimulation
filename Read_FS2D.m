% Read ASCII output
% Clear workspace, close windows, clear cli
clear all; close all; clc;
% Open file 1
FileID = fopen('Gaussian_test-00000000-0000.dat');
% Read number of cells on x axys
IMAX = fscanf(FileID,'%d',1);
% Read number of cells on y axys1
JMAX = fscanf(FileID,'%d',1);
% Read data
x   = fscanf(FileID,'%f \n',IMAX);
y   = fscanf(FileID,'%f \n',JMAX);
eta = zeros(IMAX,JMAX);
u   = zeros(IMAX,JMAX);
v   = zeros(IMAX,JMAX);
for i=1:IMAX
    eta(i,:) = fscanf(FileID,'%f \n',JMAX);
end
for i=1:IMAX
    u(i,:)   = fscanf(FileID,'%f \n',JMAX);
end
for i=1:IMAX
    v(i,:)   = fscanf(FileID,'%f \n',JMAX);
end

% Plot data
figure(1);
mesh(x,y,eta)
title('Free surface elevation')

figure(2);
mesh(x,y,u)
title('Velocity on x axys')

figure(3);
mesh(x,y,v)
title('Velocity on y axys')

% Open file 2
FileID = fopen('Gaussian_test-00000006-0000.dat');
% Read number of cells on x axys
IMAX = fscanf(FileID,'%d',1);
% Read number of cells on y axys1
JMAX = fscanf(FileID,'%d',1);
% Read data
x   = fscanf(FileID,'%f \n',IMAX);
y   = fscanf(FileID,'%f \n',JMAX);
eta = zeros(IMAX,JMAX);
u   = zeros(IMAX,JMAX);
v   = zeros(IMAX,JMAX);
for i=1:IMAX
    eta(i,:) = fscanf(FileID,'%f \n',JMAX);
end
for i=1:IMAX
    u(i,:)   = fscanf(FileID,'%f \n',JMAX);
end
for i=1:IMAX
    v(i,:)   = fscanf(FileID,'%f \n',JMAX);
end

% Plot data
figure(4);
mesh(x,y,eta)
title('Free surface elevation')

figure(5);
mesh(x,y,u)
title('Velocity on x axys')

figure(6);
mesh(x,y,v)
title('Velocity on y axys')

% Open file 3
FileID = fopen('Gaussian_test-00000012-0000.dat');
% Read number of cells on x axys
IMAX = fscanf(FileID,'%d',1);
% Read number of cells on y axys1
JMAX = fscanf(FileID,'%d',1);
% Read data
x   = fscanf(FileID,'%f \n',IMAX);
y   = fscanf(FileID,'%f \n',JMAX);
eta = zeros(IMAX,JMAX);
u   = zeros(IMAX,JMAX);
v   = zeros(IMAX,JMAX);
for i=1:IMAX
    eta(i,:) = fscanf(FileID,'%f \n',JMAX);
end
for i=1:IMAX
    u(i,:)   = fscanf(FileID,'%f \n',JMAX);
end
for i=1:IMAX
    v(i,:)   = fscanf(FileID,'%f \n',JMAX);
end

% Plot data
figure(7);
mesh(x,y,eta)
title('Free surface elevation')

figure(8);
mesh(x,y,u)
title('Velocity on x axys')

figure(9);
mesh(x,y,v)
title('Velocity on y axys')