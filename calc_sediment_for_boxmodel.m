%------------------------------------------------------------
% Calculates the sizes and surface area for the boxes of a 
% 12-box model for the global biogeochemical cycles in the ocean. 
% also calculates the area of seafloor above a depth of XXXm
% that is contained in each box, to estimate the sediment
% flux of iron into the box. 
%------------------------------------------------------------

% this is where I have the sewater equation of state for calculating sigma
addpath /Users/cvoelker/matlab/tools/whoi

% this is where I have the hydrographic data
datapath='/Users/cvoelker/data/climatologies/WOA09/';

a = read_one_netcdf([datapath,'temperature_annual_1deg.nc']);
temp = a.t_an;
fill_value = a.attributes.t_an.FillValue; 
lon = a.lon;
dlon = a.lon_bnds(2,:) - a.lon_bnds(1,:);
lat = a.lat;
dlat = a.lat_bnds(2,:) - a.lat_bnds(1,:);
depth = a.depth;

b = read_one_netcdf([datapath,'salinity_annual_1deg.nc']);
salt = b.s_an;

% calculate sigma0 for separating AAIW
[svol,sigma0] = swstate(salt, temp, zeros(size(salt)));
% calculate sigma4 for separating AABW and LCDW
[svol,sigma4] = swstate(salt, temp, 4000*ones(size(salt)));

nx = length(lon);
ny = length(lat);
nz = length(depth);

depth_bnds = [depth(1) ; (depth(1:end-1)+depth(2:end))/2 ; ...
    depth(end)*3/2 - depth(end-1)/2];
depth_int = diff(depth_bnds);

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

Atlant_mask  = ismember(basincode,[1, 30, 31, 36, 37, 48, 49, 50, 52, 58]);
Indopac_mask = ismember(basincode,[2, 3, 12, 32, 34, 35, 38, 42, 43, 44, 45, 46, 51, 56]);
SOcean_mask  = (basincode==10);

% split Southern Ocean into an Atlantic and an Indopacific part
SO_Atl_mask = SOcean_mask;
ii = find(lon>20 & lon<290);
SO_Atl_mask(ii,:,:) = 0;
SO_Pac_mask = SOcean_mask - SO_Atl_mask;

% add SO partial masks to IndoPac and Atl masks
Atlant_mask_all = Atlant_mask + SO_Atl_mask;
Indopac_mask_all = Indopac_mask + SO_Pac_mask;

%------------------------------------------------------------
% read data for global earth topography in 2 arcminutes
%------------------------------------------------------------

datapath = '/Users/cvoelker/data/climatologies/ETOPO/';
test = read_one_netcdf([datapath,'ETOPO2v2g_f4.nc']);

x_etopo = test.Longitude;
y_etopo = test.Latitude;
nxe = length(x_etopo);
nye = length(y_etopo);
rose = (-1.0)*test.z;

clear test

% now create a 3-dimensional field of the same dimension as the 
% WOA files, and that contains how much sediment area is contained
% in each box.

sed_area_3d = zeros([nx,ny,nz]);
r_earth = 6430.0e3; % mean earth radius
cosy = cosd(y_etopo);
dx = x_etopo(2) - x_etopo(1);
dy = y_etopo(2) - y_etopo(1);
dry = 2*pi*r_earth * dy / 360.0;
drx = 2*pi*r_earth * cosy * dx / 360.0;
darea = drx * dry;

for i=1:nxe,
  xb = x_etopo(i);
  if xb<0, xb = xb+360; end  % convert latitudes<0
  ib = ceil(xb);             % this is the longitude index in 3-d
  for j=2:nye,
    yb = y_etopo(j);
    jb = ceil(yb)+90;
    
    etopo_depth = rose(i,j);
    etopo_area = darea(j);
    k = sum(etopo_depth>depth_bnds);
%    sed_area_3d(ib,jb,k) = sed_area_3d(ib,jb,k) + etopo_area;
    if (k>0 && k<=nz),
      sed_area_3d(ib,jb,k) = sed_area_3d(ib,jb,k) + etopo_area;
    end
  end
end

%------------------------------------------------------------
% calculate sediment iron flux from sediment area times a maximum flux,
% times a depth dependency. The depth dependency is taken from the 
% PISCES manual, Aumont et al. 2015.
%------------------------------------------------------------

zfsed = (min(8,(depth/500).^(-1.5)));
zeta_fsed = -0.9543+0.7662*log(zfsed) -  0.235*(log(zfsed)).^2;
Fsed = min(1,exp(zeta_fsed)/0.5);
Fmax = 2.0e-6 * 360.0; % mol Fe/m^2/yr

fe_flux_3d = zeros(size(sed_area_3d));
for k=1:nz,
  fe_flux_3d(:,:,k) = sed_area_3d(:,:,k)*Fmax*Fsed(k); % mol Fe/yr
end

%------------------------------------------------------------
% to calculate the volumes of specific water masses we need
% the areas and the vertical boundaries of the boxes 
%------------------------------------------------------------

% calculate the area of an ocean box in m^2
area = zeros(nx,ny);
for i=1:nx
    for j=1:ny
        area(i,j) = (2*pi*r_earth * dlat(j)/ 360.0) * ...
            (2*pi*r_earth * cosd(lat(j)) * dlon(i) / 360.0);
    end
end

% now calculate the volumes of boxes
volumes = zeros(nx,ny,nz);
for k=1:nz
    volumes(:,:,k) = area(:,:) * depth_int(k);
end

%------------------------------------------------------------
% now start calculating the volumes of boxes of interest. For
% that first construct masks for latitude regions and for
% density boundaries
%------------------------------------------------------------

zero_mask = zeros(size(volumes));

lat_mask_north = zero_mask;
jj = find(lat >=24);
lat_mask_north(:,jj,:) = 1;

lat_mask_equat = zero_mask;
jj = find(lat<24 & lat>=-32);
lat_mask_equat(:,jj,:) = 1;

lat_mask_south = zero_mask;
jj = find(lat<-32);
lat_mask_south(:,jj,:) = 1;

lat_mask_SO_lo = zero_mask;
jj = find(lat<-60);
lat_mask_SO_lo(:,jj,:) = 1;

depth_mask = zero_mask;
kk = find(depth<=100);
depth_mask(:,:,kk) = 1;

sigma0_mask = zero_mask;
sigma0_mask(sigma0<27.3) = 1;

sigma4_mask = zero_mask;
sigma4_mask(sigma4>45.86) = 1;

%------------------------------------------------------------
% Now construct the boxes as combination of masks
%------------------------------------------------------------

% surface boxes
SO_surf = depth_mask .* lat_mask_south .* ...
    (Atlant_mask_all + Indopac_mask_all);
EqPac_surf = depth_mask .* lat_mask_equat .* Indopac_mask;
NPac_surf = depth_mask .* lat_mask_north .* Indopac_mask;
EqAtl_surf = depth_mask .* lat_mask_equat .* Atlant_mask;
NAtl_surf = depth_mask .* lat_mask_north .* Atlant_mask;

% deep SO box
SO_deep = ~depth_mask .* lat_mask_SO_lo .* ...
    (Atlant_mask_all + Indopac_mask_all);

% AAIW boxes
AAIW_Pac = Indopac_mask_all .* sigma0_mask .* (~depth_mask) .* ...
    (~lat_mask_SO_lo);
AAIW_Atl = Atlant_mask_all  .* sigma0_mask .* (~depth_mask) .* ...
    (~lat_mask_SO_lo);

% deep boxes
Atl_deep = Atlant_mask_all .* (~lat_mask_SO_lo) .* ...
    ~depth_mask .* ~sigma0_mask .* ~sigma4_mask;
Atl_bottom = Atlant_mask_all .* (~lat_mask_SO_lo) .* ...
    ~depth_mask .* sigma4_mask;

Pac_deep = Indopac_mask_all .* (~lat_mask_SO_lo) .* ...
    ~depth_mask .* ~sigma0_mask .* ~sigma4_mask;
Pac_bottom = Indopac_mask_all .* (~lat_mask_SO_lo) .* ...
    ~depth_mask .* sigma4_mask;

%------------------------------------------------------------
% calculate the volumes of the boxes and
% the surface area of the surface boxes
%------------------------------------------------------------

boxvolumes = zeros(12,1);
boxvolumes(1) = sum(sum(sum( volumes .* NPac_surf)));
boxvolumes(2) = sum(sum(sum( volumes .* EqPac_surf)));
boxvolumes(3) = sum(sum(sum( volumes .* SO_surf)));
boxvolumes(4) = sum(sum(sum( volumes .* EqAtl_surf)));
boxvolumes(5) = sum(sum(sum( volumes .* NAtl_surf)));
boxvolumes(6) = sum(sum(sum( volumes .* AAIW_Pac)));
boxvolumes(7) = sum(sum(sum( volumes .* AAIW_Atl)));
boxvolumes(8) = sum(sum(sum( volumes .* Pac_deep)));
boxvolumes(9) = sum(sum(sum( volumes .* SO_deep)));
boxvolumes(10) = sum(sum(sum( volumes .* Atl_deep)));
boxvolumes(11) = sum(sum(sum( volumes .* Pac_bottom)));
boxvolumes(12) = sum(sum(sum( volumes .* Atl_bottom)));

boxareas = zeros(5,1);
boxareas(1) = sum( sum( area .* NPac_surf(:,:,1)));
boxareas(2) = sum( sum( area .* EqPac_surf(:,:,1)));
boxareas(3) = sum( sum( area .* SO_surf(:,:,1)));
boxareas(4) = sum( sum( area .* EqAtl_surf(:,:,1)));
boxareas(5) = sum( sum( area .* NAtl_surf(:,:,1)));

%------------------------------------------------------------
% now sum up the total Fe source in each of the boxes 
%------------------------------------------------------------

box_fe = zeros(12,1);
box_fe(1) = sum(sum(sum( fe_flux_3d  .* NPac_surf)));
box_fe(2) = sum(sum(sum( fe_flux_3d  .* EqPac_surf)));
box_fe(3) = sum(sum(sum( fe_flux_3d  .* SO_surf)));
box_fe(4) = sum(sum(sum( fe_flux_3d  .* EqAtl_surf)));
box_fe(5) = sum(sum(sum( fe_flux_3d  .* NAtl_surf)));
box_fe(6) = sum(sum(sum( fe_flux_3d  .* AAIW_Pac)));
box_fe(7) = sum(sum(sum( fe_flux_3d  .* AAIW_Atl)));
box_fe(8) = sum(sum(sum( fe_flux_3d  .* Pac_deep)));
box_fe(9) = sum(sum(sum( fe_flux_3d  .* SO_deep)));
box_fe(10) = sum(sum(sum( fe_flux_3d .* Atl_deep)));
box_fe(11) = sum(sum(sum( fe_flux_3d .* Pac_bottom)));
box_fe(12) = sum(sum(sum( fe_flux_3d .* Atl_bottom)));

boxfe_change_sed = box_fe ./ boxvolumes * 1.0e6; % (from mol/m^3 to nmol/L)

fprintf('total Fe flux (mol/yr) into the boxes and induced rate of change (nmol/L/yr)\n')
for k=1:12,
    fprintf('%10.2e %10.2e\n', box_fe(k), boxfe_change_sed(k))
end

return



