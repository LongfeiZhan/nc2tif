clear
clc
% 2010-2013

% information of nc
info = ncinfo('1.nc');
% read variables
lon = ncread('1.nc','longitude');
lat = ncread('1.nc','latitude');
[mx,my] = meshgrid(lon,lat); % meshing the longitude and latitude

lev = ncread('1.nc','level');
time = ncread('1.nc','time');
% exchange date to a easier formate
t0 = datetime(1900,1,1);
date_yyymmdd = t0 + double(time(:))/24;
t_str = datestr(date_yyymmdd,'yyyymmdd');   %datetime??????char??

% Get time conversion lines
r = 1;
for i = 1 : length(t_str)-1
    if strcmp(t_str(i,1:4),t_str(i+1,1:4))
        continue
    else
        mark(r,1) = i;
        r = r + 1;
    end
end
mark(r,1) = length(t_str);
mark = [0;mark]; % time conversion lines

%% 2010
tem2010 = ncread('1.nc','t',[1,1,1,mark(1)+1],[Inf,Inf,Inf,mark(2)]);
tem2010_ave = mean(tem2010,4);

tem_ave_975 = flipud(tem2010_ave(:,:,1)');  % 975hpa tem ave
rasterSize = size(tem_ave_975);        % size
R = georasterref('RasterSize', rasterSize,'Latlim', [double(min(lat)),double(max(lat))],'Lonlim',[double(min(lon)),double(max(lon))]);
geotiffwrite('2010up.tif',tem_ave_975,R);  % writting .tif file

tem_ave_1000 = flipud(tem2010_ave(:,:,2)');  % 1000hpa tem ave
rasterSize = size(tem_ave_1000);        % size
R = georasterref('RasterSize', rasterSize,'Latlim', [double(min(lat)),double(max(lat))],'Lonlim',[double(min(lon)),double(max(lon))]);
geotiffwrite('2010down.tif',tem_ave_1000,R); % writting .tif file

tem2010_inversion = tem2010_ave(:,:,1) - tem2010_ave(:,:,2); % inversion
for i = 1 : size(tem2010_inversion,1)
    for j = 1 : size(tem2010_inversion,2)
        if tem2010_inversion(i,j) < 0
            tem2010_inversion(i,j) = 0; % if inversion < 0,then inversion = 0
        else
            disp(['2010','?',num2str(i),'??',num2str(j),'?????'])
        end
    end
end
tem2010_inv = flipud(tem2010_inversion');  % 1000hpa tem ave
rasterSize = size(tem2010_inv);        % size
R = georasterref('RasterSize', rasterSize,'Latlim', [double(min(lat)),double(max(lat))],'Lonlim',[double(min(lon)),double(max(lon))]);
geotiffwrite('2010inversion.tif',tem2010_inv,R);  % writting .tif file

clear tem2010  % Release memory

%% 2011
tem2011 = ncread('1.nc','t',[1,1,1,mark(2)+1],[Inf,Inf,Inf,mark(3)-mark(2)]);
tem2011_ave = mean(tem2011,4);
tem_ave_975 = flipud(tem2011_ave(:,:,1)');  % 975hpa tem ave
rasterSize = size(tem_ave_975);        % size
R = georasterref('RasterSize', rasterSize,'Latlim', [double(min(lat)),double(max(lat))],'Lonlim',[double(min(lon)),double(max(lon))]);
geotiffwrite('2011up.tif',tem_ave_975,R);  % writting .tif file

tem_ave_1000 = flipud(tem2011_ave(:,:,2)');  % 1000hpa tem ave
rasterSize = size(tem_ave_1000);        % size
R = georasterref('RasterSize', rasterSize,'Latlim', [double(min(lat)),double(max(lat))],'Lonlim',[double(min(lon)),double(max(lon))]);
geotiffwrite('2011down.tif',tem_ave_1000,R); % writting .tif file

tem2011_inversion = tem2011_ave(:,:,1) - tem2011_ave(:,:,2); % inversion
for i = 1 : size(tem2011_inversion,1)
    for j = 1 : size(tem2011_inversion,2)
        if tem2011_inversion(i,j) < 0
            tem2011_inversion(i,j) = 0; % if inversion < 0,then inversion = 0
        else
            disp(['2011','?',num2str(i),'??',num2str(j),'?????'])
        end
    end
end
tem2011_inv = flipud(tem2011_inversion');  % 1000hpa tem ave
rasterSize = size(tem2011_inv);        % size
R = georasterref('RasterSize', rasterSize,'Latlim', [double(min(lat)),double(max(lat))],'Lonlim',[double(min(lon)),double(max(lon))]);
geotiffwrite('2011inversion.tif',tem2011_inv,R);  % writting .tif file

clear tem2011 % Release memory

%% 2012
tem2012 = ncread('1.nc','t',[1,1,1,mark(3)+1],[Inf,Inf,Inf,mark(4)-mark(3)]);
tem2012_ave = mean(tem2012,4);

tem_ave_975 = flipud(tem2012_ave(:,:,1)');  % 975hpa tem ave
rasterSize = size(tem_ave_975);        % size
R = georasterref('RasterSize', rasterSize,'Latlim', [double(min(lat)),double(max(lat))],'Lonlim',[double(min(lon)),double(max(lon))]);
geotiffwrite('2012up.tif',tem_ave_975,R);  % writting .tif file

tem_ave_1000 = flipud(tem2012_ave(:,:,2)');  % 1000hpa tem ave
rasterSize = size(tem_ave_1000);        % size
R = georasterref('RasterSize', rasterSize,'Latlim', [double(min(lat)),double(max(lat))],'Lonlim',[double(min(lon)),double(max(lon))]);
geotiffwrite('2012down.tif',tem_ave_1000,R); % writting .tif file

tem2012_inversion = tem2012_ave(:,:,1) - tem2012_ave(:,:,2); % inversion
for i = 1 : size(tem2012_inversion,1)
    for j = 1 : size(tem2012_inversion,2)
        if tem2012_inversion(i,j) < 0
            tem2012_inversion(i,j) = 0; % if inversion < 0,then inversion = 0
        else
            disp(['2012','?',num2str(i),'??',num2str(j),'?????'])
        end
    end
end
tem2012_inv = flipud(tem2012_inversion');  % 1000hpa tem ave
rasterSize = size(tem2012_inv);        % size
R = georasterref('RasterSize', rasterSize,'Latlim', [double(min(lat)),double(max(lat))],'Lonlim',[double(min(lon)),double(max(lon))]);
geotiffwrite('2012inversion.tif',tem2012_inv,R);  % writting .tif file

clear tem2012

%% 2013
tem2013 = ncread('1.nc','t',[1,1,1,mark(4)+1],[Inf,Inf,Inf,mark(5)-mark(4)]);
tem2013_ave = mean(tem2013,4);

tem_ave_975 = flipud(tem2013_ave(:,:,1)');  % 975hpa tem ave
rasterSize = size(tem_ave_975);        % size
R = georasterref('RasterSize', rasterSize,'Latlim', [double(min(lat)),double(max(lat))],'Lonlim',[double(min(lon)),double(max(lon))]);
geotiffwrite('2013up.tif',tem_ave_975,R);  % writting .tif file

tem_ave_1000 = flipud(tem2013_ave(:,:,2)');  % 1000hpa tem ave
rasterSize = size(tem_ave_1000);        % size
R = georasterref('RasterSize', rasterSize,'Latlim', [double(min(lat)),double(max(lat))],'Lonlim',[double(min(lon)),double(max(lon))]);
geotiffwrite('2013down.tif',tem_ave_1000,R); % writting .tif file

tem2013_inversion = tem2013_ave(:,:,1) - tem2013_ave(:,:,2); % inversion
for i = 1 : size(tem2013_inversion,1)
    for j = 1 : size(tem2013_inversion,2)
        if tem2013_inversion(i,j) < 0
            tem2013_inversion(i,j) = 0; % if inversion < 0,then inversion = 0
        else
            disp(['2013','?',num2str(i),'??',num2str(j),'?????'])
        end
    end
end
tem2013_inv = flipud(tem2013_inversion');  % 1000hpa tem ave
rasterSize = size(tem2013_inv);        % size
R = georasterref('RasterSize', rasterSize,'Latlim', [double(min(lat)),double(max(lat))],'Lonlim',[double(min(lon)),double(max(lon))]);
geotiffwrite('2013inversion.tif',tem2013_inv,R);  % writting .tif file

clear tem2013