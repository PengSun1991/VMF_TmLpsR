load('orography.mat'); %load the orography.mat file, which contains the height of each VMF1 grid point
load('btmodel.mat');%load the tm lapse-rate model
lat=44.5/180*pi; % latitude in radian
lon=116.4/180*pi; %longitude in radian
h_ell=100; %height in m
mjd=5.848484409722220e+04;%2019 1 1 20 15 30
indir_VMF1_grid='.\STD_FC';
VMF1_grid_file=[];
[tm,VMF1_grid_file] = vmf1_grid_tm (indir_VMF1_grid,orography,VMF1_grid_file,mjd,lat,lon,h_ell,btmodel);