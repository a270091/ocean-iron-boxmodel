%------------------------------------------------------------
% Calculates the sizes and surface area for the boxes of a 
% 12-box model for the global biogeochemical cycles in the ocean. 
% also calculates the average values of phosphate and the total 
% of export production 
%------------------------------------------------------------

% this is where I have the routine 'read_one_netcdf.m'
addpath /Users/cvoelker/matlab/tools/MITgcm

% this is where I have the sewater equation of state for calculating sigma
addpath /Users/cvoelker/matlab/tools/whoi

% this is where I have the data
datapath='/Users/cvoelker/data/climatologies/WOA09/';

a = read_one_netcdf([datapath,'phosphate_annual_1deg.nc']);
po4 = a.p_an;
fill_value = a.attributes.p_an.FillValue; 
lon = a.lon;
dlon = a.lon_bnds(2,:) - a.lon_bnds(1,:);
lat = a.lat;
dlat = a.lat_bnds(2,:) - a.lat_bnds(1,:);
depth = a.depth;

b = read_one_netcdf([datapath,'silicate_annual_1deg.nc']);
sil = b.i_an;
b = read_one_netcdf([datapath,'temperature_annual_1deg.nc']);
temp = b.t_an;
b = read_one_netcdf([datapath,'salinity_annual_1deg.nc']);
salt = b.s_an;

% calculate sigma0 for separating AAIW
[svol,sigma0] = swstate(salt, temp, zeros(size(salt)));
% calculate sigma4 for separating AABW and LCDW
[svol,sigma4] = swstate(salt, temp, 4000*ones(size(salt)));

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
% read export production from Schlitzer; attention, the export
% map starts at 20.5E, not at 0.5E. So we need to shift things a bit
%------------------------------------------------------------

datapath='/Users/cvoelker/data/climatologies/export/';
c = read_one_netcdf([datapath,'export_schlitzer.nc']);
export = c.POC_SCHLITZER;
export(export==c.attributes.POC_SCHLITZER.missing_value) = 0;

% do the shift
jj = find(c.ETOPO60X>360);
jj2 = find(c.ETOPO60X<360);
export2 = [export(jj,:); export(jj2,:)];

%------------------------------------------------------------
% to calculate the volumes of specific water masses we need
% the areas and the vertical boundaries of the boxes 
%------------------------------------------------------------

nx = length(lon);
ny = length(lat);
nz = length(depth);

depth_bnds = [depth(1) ; (depth(1:end-1)+depth(2:end))/2 ; ...
    depth(end)*3/2 - depth(end-1)/2];
depth_int = diff(depth_bnds);

% calculate the area of an ocean box in m^2
r_earth = 6430.0e3; % mean earth radius
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
% calculate the volumes, and the average PO4 concentration of the boxes and
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

for k=1:12,
    fprintf('Box %i, Volume %8.3e m**3\n',k,boxvolumes(k));
end
fprintf('\n');

boxareas = zeros(5,1);
boxareas(1) = sum( sum( area .* NPac_surf(:,:,1)));
boxareas(2) = sum( sum( area .* EqPac_surf(:,:,1)));
boxareas(3) = sum( sum( area .* SO_surf(:,:,1)));
boxareas(4) = sum( sum( area .* EqAtl_surf(:,:,1)));
boxareas(5) = sum( sum( area .* NAtl_surf(:,:,1)));

for k=1:5,
    fprintf('Box %i, Area %8.3e m**2\n',k,boxareas(k));
end
fprintf('\n');

% now average PO4, Silicate, Temperature, and Salinity over the boxes

boxpo4 = zeros(12,1);
boxpo4(1) = sum(sum(sum( po4 .* volumes .* NPac_surf)));
boxpo4(2) = sum(sum(sum( po4 .* volumes .* EqPac_surf)));
boxpo4(3) = sum(sum(sum( po4 .* volumes .* SO_surf)));
boxpo4(4) = sum(sum(sum( po4 .* volumes .* EqAtl_surf)));
boxpo4(5) = sum(sum(sum( po4 .* volumes .* NAtl_surf)));
boxpo4(6) = sum(sum(sum( po4 .* volumes .* AAIW_Pac)));
boxpo4(7) = sum(sum(sum( po4 .* volumes .* AAIW_Atl)));
boxpo4(8) = sum(sum(sum( po4 .* volumes .* Pac_deep)));
boxpo4(9) = sum(sum(sum( po4 .* volumes .* SO_deep)));
boxpo4(10) = sum(sum(sum( po4 .* volumes .* Atl_deep)));
boxpo4(11) = sum(sum(sum( po4 .* volumes .* Pac_bottom)));
boxpo4(12) = sum(sum(sum( po4 .* volumes .* Atl_bottom)));
boxpo4 = boxpo4 ./ boxvolumes;

for k=1:12,
    fprintf('%4.2f\n', boxpo4(k))
end

exports = zeros(5,1);
exports(1) = sum( sum( export2 .* area .* NPac_surf(:,:,1)));
exports(2) = sum( sum( export2 .* area .* EqPac_surf(:,:,1)));
exports(3) = sum( sum( export2 .* area .* SO_surf(:,:,1)));
exports(4) = sum( sum( export2 .* area .* EqAtl_surf(:,:,1)));
exports(5) = sum( sum( export2 .* area .* NAtl_surf(:,:,1)));
exports_Pg = exports .* 1.0e-15;
exports_per_area = exports ./ boxareas / 12;

fprintf('Export (PgC/yr)\n')
for k=1:5,
    fprintf('%4.2f\n', exports_Pg(k))
end
fprintf('C flux (mol C/(m^2 yr))\n')
for k=1:5,
    fprintf('%4.2f\n', exports_per_area(k))
end


return

