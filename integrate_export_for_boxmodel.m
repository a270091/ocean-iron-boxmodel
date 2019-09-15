clear all
%-----------------------------------------------
% add path to read netCDF files
%-----------------------------------------------
ismac = strcmp(computer('arch'),'maci64');
if (ismac==1),
  addpath /Users/cvoelker/matlab/tools/MITgcm/
  addpath /Users/cvoelker/matlab/tools/m_map/
end

%------------------------------------------------------------
% read export production from Schlitzer or Laws; attention, the export
% map starts at 20.5E, not at 0.5E. So we need to shift things a bit
%------------------------------------------------------------

% which export? 
iex = 2;

datapath='/Users/cvoelker/data/climatologies/export/';
if (iex==1),
  c = read_one_netcdf([datapath,'export_schlitzer.nc']);
  export = c.POC_SCHLITZER;
  export(export==c.attributes.POC_SCHLITZER.missing_value) = 0;
elseif(iex==2),
  c = read_one_netcdf([datapath,'laws_export.nc']);
  export = c.EXP;
  export(export==c.attributes.EXP.missing_value) = 0;
else
  c = read_one_netcdf([datapath,'eppley_export.nc']);
  export = c.EXP;
  export(export==c.attributes.EXP.missing_value) = 0;
end

% do the shift
jj = find(c.ETOPO60X>360);
jj2 = find(c.ETOPO60X<360);
export2 = [export(jj,:); export(jj2,:)];

lon = c.ETOPO60X;
dlon = lon(2) - lon(1);
lat = c.ETOPO60Y;
dlat = lat(2) - lat(1);
% depth = a.depth;

%------------------------------------------------------------
% to calculate the volumes of specific water masses we need
% the areas and the vertical boundaries of the boxes 
%------------------------------------------------------------

nx = length(lon);
ny = length(lat);
%nz = length(depth);

%depth_bnds = [depth(1) ; (depth(1:end-1)+depth(2:end))/2 ; ...
%    depth(end)*3/2 - depth(end-1)/2];
%depth_int = diff(depth_bnds);

% calculate the area of an ocean box in m^2
r_earth = 6430.0e3; % mean earth radius
area = zeros(nx,ny);
for i=1:nx
    for j=1:ny
        area(i,j) = (2*pi*r_earth * dlat/ 360.0) * ...
            (2*pi*r_earth * cosd(lat(j)) * dlon / 360.0);
    end
end

% now calculate the volumes of boxes
%volumes = zeros(nx,ny,nz);
%for k=1:nz
%    volumes(:,:,k) = area(:,:) * depth_int(k);
%end

%-----------------------------------------------------------
% read the ocean WOA basin code data
% basincode==1:  Atlantic
% basincode==2:  Pacific
% basincode==3:  Indian Ocean
% basincode==10: Southern Ocean
%-----------------------------------------------------------

datapath='/Users/cvoelker/data/climatologies/export/';
c = read_one_netcdf([datapath,'oceanbasins.nc']);
basincode = c.basin;
lon2 = c.X;
lat2 = c.Y;

Atlant_mask  = ismember(basincode,[1, 30, 31, 36, 37, 48, 49, 50, 52, 58]);
Indopac_mask = ismember(basincode,[2, 3, 12, 32, 34, 35, 38, 42, 43, 44, 45, 46, 51, 56]);
SOcean_mask  = (basincode==10);

% split Southern Ocean into an Atlantic and an Indopacific part
SO_Atl_mask = SOcean_mask;
ii = find(lon2>20 & lon2<290);
SO_Atl_mask(ii,:,:) = 0;
SO_Pac_mask = SOcean_mask - SO_Atl_mask;

% add SO partial masks to IndoPac and Atl masks
Atlant_mask_all = Atlant_mask + SO_Atl_mask;
Indopac_mask_all = Indopac_mask + SO_Pac_mask;

Atlant_mask_surf = squeeze(Atlant_mask_all(:,:,1));
Indopac_mask_surf = squeeze(Indopac_mask_all(:,:,1));

zero_mask = zeros(size(Atlant_mask_surf));

lat_mask_north = zero_mask;
jj = find(lat2 >=24);
lat_mask_north(:,jj) = 1;

lat_mask_equat = zero_mask;
jj = find(lat2<24 & lat2>=-32);
lat_mask_equat(:,jj) = 1;

lat_mask_south = zero_mask;
jj = find(lat2<-32);
lat_mask_south(:,jj) = 1;

lat_mask_deepsouth = zero_mask;
jj = find(lat<-60);
lat_mask_deepsouth(:,jj,:) = 1;

box1_mask = (Indopac_mask_surf .* lat_mask_north);
box2_mask = (Indopac_mask_surf .* lat_mask_equat);
box3_mask = ((Indopac_mask_surf + Atlant_mask_surf) .* lat_mask_south);
box3_mask_s = ((Indopac_mask_surf + Atlant_mask_surf) .* lat_mask_deepsouth);
box3_mask_a = (Atlant_mask_surf .* (lat_mask_south - lat_mask_deepsouth));
box3_mask_p = (Indopac_mask_surf .* (lat_mask_south - lat_mask_deepsouth));
box4_mask = (Atlant_mask_surf .* lat_mask_equat);
box5_mask = (Atlant_mask_surf .* lat_mask_north);

exports = zeros(5,1);
boxareas = zeros(5,1);
exports(1) = sum( sum( export2 .* area .* box1_mask));
exports(2) = sum( sum( export2 .* area .* box2_mask));
exports(3) = sum( sum( export2 .* area .* box3_mask));
exports(4) = sum( sum( export2 .* area .* box4_mask));
exports(5) = sum( sum( export2 .* area .* box5_mask));
exp3_s = sum( sum( export2 .* area .* box3_mask_s));
exp3_a = sum( sum( export2 .* area .* box3_mask_a));
exp3_p = sum( sum( export2 .* area .* box3_mask_p));
boxareas(1) = sum( sum( area .* box1_mask));
boxareas(2) = sum( sum( area .* box2_mask));
boxareas(3) = sum( sum( area .* box3_mask));
boxareas(4) = sum( sum( area .* box4_mask));
boxareas(5) = sum( sum( area .* box5_mask));

exports_Pg = exports .* 1.0e-15;
exports_per_area = exports ./ boxareas / 12;

fprintf('Export (PgC/yr)\n')
for k=1:5,
  fprintf('%4.2f\n', exports_Pg(k))
  if k==3,
    x_s = exp3_s/exports(3);
    x_a = exp3_a/exports(3);
    x_p = exp3_p/exports(3);
    fprintf('split into 3 parts: S %4.2f A %4.2f P %4.2f \n',x_s,x_a,x_p)
  end
end
fprintf('C flux (mol C/(m^2 yr))\n')
for k=1:5,
    fprintf('%4.2f\n', exports_per_area(k))
end


return
