function [tm,VMF1_grid_file] = vmf1_grid_tm (indir_VMF1_grid,orography_ell,VMF1_grid_file,mjd,lat,lon,h_ell,btmodel)
%
% vmf1_grid_tm.m
% This routine determines weighted mean temperature from the gridded Compact VMF1 files, as available from:
% Forecast:
%   http://vmf.geo.tuwien.ac.at/trop_products/GRID/2.5x2/VMF1/STD_FC
% Operational:
%   http://vmf.geo.tuwien.ac.at/trop_products/GRID/2.5x2/VMF1/STD_OP
%
% This file is modified from the official function 'vmf1_grid.m' :
%   https://vmf.geo.tuwien.ac.at/codes/vmf1_grid.m
% The authors greatly appreciate the Tu Wien for providing the code that written by Daniel Landskron.

% INPUT:
%         o indir_VMF1_grid ... input directory where the yearly subdivided VMF1 gridded files are stored
%         o orography_ell ...   ellipsoidal height of the grid points
%         o VMF1_grid_file: ... cell containing filenames, VMF1 data and the orography, which is always passed with the function, must be set to '[]' by the user in the initial run
%         o mjd ............... modified Julian date
%         o lat ............... ellipsoidal latitude in radians
%         o lon ............... ellipsoidal longitude in radians
%         o h_ell ............. ellipsoidal height in meters
%         o btmodel ............the Tm lapse-rate model 
%
% OUTPUT:
%         o tm  ............... tm, valid at h_ell
%         o VMF1_grid_file: ... cell containing filenames, VMF1 data and the orography, which is always passed with the function, must be set to '[]' by the user in the initial run
%
% -------------------------------------------------------------------------
%
% written by Peng Sun (2024-07-17)
% peng_sun@cumt.edu.cn
%
% ===================================================



% save lat and lon also in degrees
lat_deg = lat*180/pi;
lon_deg = lon*180/pi;

% due to numerical issues, it might happen that the above conversion does not give exact results, e.g. in case of rad2deg(deg2rad(60)); in order to prevent this, lat_deg and lon_deg are rounded to the 10th decimal place
lat_deg = round(lat_deg,10);
lon_deg = round(lon_deg,10);


%% (1) convert the mjd to year, month, day in order to find the correct files


% find the two surrounding epochs
if mod(mjd,0.25)==0
    mjd_all = mjd;
else
    mjd_int = floor(mjd*4)/4 : 0.25 : ceil(mjd*4)/4;
    mjd_all = [mjd mjd_int];
end


hour = floor((mjd_all-floor(mjd_all))*24);   % get hours
minu = floor((((mjd_all-floor(mjd_all))*24)-hour)*60);   % get minutes
sec = (((((mjd_all-floor(mjd_all))*24)-hour)*60)-minu)*60;   % get seconds

% change secs, min hour whose sec==60
minu(sec==60) = minu(sec==60)+1;
hour(minu==60) = hour(minu==60)+1;
mjd_all(hour==24)=mjd_all(hour==24)+1;

% calc jd (yet wrong for hour==24)
jd_all = mjd_all+2400000.5;

% integer Julian date
jd_all_int = floor(jd_all+0.5);

aa = jd_all_int+32044;
bb = floor((4*aa+3)/146097);
cc = aa-floor((bb*146097)/4);
dd = floor((4*cc+3)/1461);
ee = cc-floor((1461*dd)/4);
mm = floor((5*ee+2)/153);

day = ee-floor((153*mm+2)/5)+1;

month = mm+3-12*floor(mm/10);
year = bb*100+dd-4800+floor(mm/10);

epoch = (mjd_all-floor(mjd_all))*24;

% derive related VMFG filename(s)
if length(mjd_all)==1   % if the observation epoch coincides with an NWM epoch
    doy=UTC2DOY(year(1),month(1),day(1));
    filename = ['tm' num2str(year(1)-2000) sprintf('%03s',num2str(doy)) '.h' sprintf('%02s',num2str(epoch(1)))];
else
    for i_mjd = 2:length(mjd_all)
        doy=UTC2DOY(year(i_mjd),month(i_mjd),day(i_mjd));
        filename(i_mjd-1,:) = ['tm' num2str(year(i_mjd)-2000) sprintf('%03s',num2str(doy)) '.h' sprintf('%02s',num2str(epoch(i_mjd)))];
    end
end

% only positive longitude in degrees
if lon_deg < 0
    lon = lon + 2*pi;
    lon_deg = (lon_deg + 360);
end

%% (2) check if new files have to be loaded or if the overtaken ones are sufficient
if isempty(VMF1_grid_file)   % in the first run, 'VMF1_file' is always empty and the orography_ell file has to loaded
    load_new = 1;
    VMF1_grid_file{1} = filename;   % replace the empty cell by the current filenames
    VMF1_grid_file{3} = orography_ell;
    VMF1_grid_file{4} = [lat lon];
elseif strcmpi(VMF1_grid_file{1},filename)   &&   (lat > VMF1_grid_file{4}(1)   ||   (lat == VMF1_grid_file{4}(1) && lon <= VMF1_grid_file{4}(2)))   % if the current filenames are the same as in the forwarded files, and the coordinates have smaler indices as the saved ones (because the grid is only read up to the necessary indices)
    load_new = 0;
    VMF1_data_all = VMF1_grid_file{2};
else   % if new files are required, then everything must be loaded anew
    load_new = 1;
    VMF1_grid_file{1} = filename;   % replace the empty cell by the current filenames
    orography_ell = VMF1_grid_file{3};
    VMF1_grid_file{4} = [lat lon];
end
%% (3) find the indices of the 4 surrounding grid points

% get all coordinates of the grid
lat_all = 90:-2:-90;
lon_all = 0:2.5:357.5;

% find the 2 closest latitudes
lat_temp = lat_deg-lat_all;
[~,ind_lat_int(1)] = min(abs(lat_temp));
ind_lat_int(2) = ind_lat_int(1)-sign(lat_temp(ind_lat_int(1)));

% find the 2 closest longitudes
lon_temp = lon_deg-lon_all;
[~,ind_lon_int(1)] = min(abs(lon_temp));
ind_lon_int(2) = ind_lon_int(1)+sign(lon_temp(ind_lon_int(1)));

% correct indices out of range
for i_ind = 1:2
    if ind_lat_int(i_ind)>length(lat_all); ind_lat_int(i_ind) = length(lat_all);                    end
    if ind_lat_int(i_ind)<1;               ind_lat_int(i_ind) = 1;                                  end
    if ind_lon_int(i_ind)>length(lon_all); ind_lon_int(i_ind) = ind_lon_int(i_ind)-length(lon_all); end
    if ind_lon_int(i_ind)<1;               ind_lon_int(i_ind) = ind_lon_int(i_ind)+length(lon_all); end
end

% define the indices
index(1) = (ind_lat_int(1)-1)*length(lon_all)+ind_lon_int(1);
index(2) = (ind_lat_int(1)-1)*length(lon_all)+ind_lon_int(2);
index(3) = (ind_lat_int(2)-1)*length(lon_all)+ind_lon_int(1);
index(4) = (ind_lat_int(2)-1)*length(lon_all)+ind_lon_int(2);
B=90:-2:-90;
L=0:2.5:357.5;
lat_Grid=zeros(91*144,1);
lon_Grid=zeros(91*144,1);
for i=1:91
    for j=1:144
        lat_Grid((i-1)*144+j,1)=B(i);
        lon_Grid((i-1)*144+j,1)=L(j);
    end
end

%% (4) read the correct data and perform a linear time interpolation from the surrounding two epochs
% read in with textscan, but only up to maximum index, everything before will be treated as headerlines
if load_new == 1

    for i_file = 1:size(filename,1)

        % read the files and collect the data
        if length(mjd_all)==1   % if the observation epoch coincides with an NWM epoch
            path=[indir_VMF1_grid '\' num2str(year(1)) '\' filename(1,:)];
        else
            path=[indir_VMF1_grid '\' num2str(year(i_file+1)) '\' filename(i_file,:)];
        end
        data=TmGridReader(path);
        data=[lat_Grid lon_Grid data'];
        VMF1_data_all{i_file}= data;
        VMF1_grid_file{2} = VMF1_data_all;   % write the VMF1 data to the forwarded variable
        VMF1_data{i_file} = VMF1_data_all{i_file}(index,:);   % reduce to the indices of the surrounding grid points

    end
else
    for i_file = 1:size(filename,1)
        VMF1_data{i_file} =VMF1_data_all{i_file}(index,:);  % reduce to the indices of the surrounding grid points
    end
end


% initialize
VMF1_data_int_h0 = zeros(4,3);
VMF1_data_int_h1 = zeros(4,3);

% do the linear time interpolation for each argument; the results are the VMF1 values for the surrounding grid points at the time of the measurement
iv_ind = 1:4;
if length(mjd_all)==1   % if the observation epoch coincides with an NWM epoch
    VMF1_data_int_h0(iv_ind,1:3) = VMF1_data{1}(iv_ind,1:3);
else   % else perform the linear interpolation
    iv_line = 1:3;
    VMF1_data_int_h0(iv_ind,iv_line) = VMF1_data{1}(iv_ind,iv_line) + (VMF1_data{2}(iv_ind,iv_line)-VMF1_data{1}(iv_ind,iv_line))*(mjd-mjd_int(1))/(mjd_int(2)-mjd_int(1));   % the appendix 'h0' means that the values are valid at zero height
end

% the first 2 columns are equal
VMF1_data_int_h1(:,1:2) = VMF1_data_int_h0(:,1:2);


%% (5) reduce Tm of the surrounding grid points to the respective height of the location
% calculate the Tm lapse rate
doy=UTC2DOY(year(1),month(1),day(1));
cosy=cos((doy-btmodel(index,3))*2*pi/365.25);
coshy=cos((doy-btmodel(index,5))*4*pi/365.25);
bt1=btmodel(index,1)+btmodel(index,2).*cosy+btmodel(index,4).*coshy;%in K/km
% perform the vertical reduction
VMF1_data_int_h1(iv_ind,3) = VMF1_data_int_h0(iv_ind,3)- bt1.*(h_ell-orography_ell(index))/1000;


%% (6) perform the bilinear interpolation
if length(unique(index)) == 1   % if the point is directly on a grid point
    tm = VMF1_data_int_h1(1,3);
else
    % bilinear interpolation (interpreted as two 1D linear interpolations for lat and lon, but programmed without subfunctions)
    % (a) linear interpolation for longitude
    if ~isequal(VMF1_data_int_h1(1,2), VMF1_data_int_h1(2,2))   % if longitude must be interpolated (that is, the point does not have a longitude on the interval [0:2.5:357.5])
        tm_lon1 = VMF1_data_int_h1(1,3) + (VMF1_data_int_h1(2,3)-VMF1_data_int_h1(1,3))*(lon_deg-VMF1_data_int_h1(1,2))/(VMF1_data_int_h1(2,2)-VMF1_data_int_h1(1,2));
        tm_lon2 = VMF1_data_int_h1(3,3) + (VMF1_data_int_h1(4,3)-VMF1_data_int_h1(3,3))*(lon_deg-VMF1_data_int_h1(3,2))/(VMF1_data_int_h1(4,2)-VMF1_data_int_h1(3,2));

    else   % if the station coincides with the longitude of the grid
        tm_lon1 = VMF1_data_int_h1(1,3);
        tm_lon2 = VMF1_data_int_h1(3,3);
    end

    % linear interpolation for latitude
    if ~isequal(VMF1_data_int_h1(1,1), VMF1_data_int_h1(3,1))
        tm = tm_lon1 + (tm_lon2-tm_lon1)*(lat_deg-VMF1_data_int_h1(1,1))/(VMF1_data_int_h1(3,1)-VMF1_data_int_h1(1,1));
    else   % if the station coincides with the latitude of the grid
        tm = tm_lon1;
    end

end


