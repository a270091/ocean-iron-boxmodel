clear all

%-----------------------------------------------
% load monthly data
% this one has an issue in 2010:
dustpath = '/Users/cvoelker/data/dust_mahowald/';
test = read_one_netcdf([dustpath,'totaldep.monthly.79-10.nc']);
% this one has a fix for 2010:
% ncload('dep.monthlywithall.79-10.nc');
dust = test.DEPJ;
lon = test.lon;
lat = test.lat;

%-------------------------
% create montly climatology
%-------------------------

[nx,ny,nt] = size(dust);
nyear = (nt / 12) - 1; % we ignore the last year (2010)

for k=1:12,
  krange = ((1:nyear)-1)*12+k;
  dust_clim_new(:,:,k)=mean(dust(:,:,krange),3);
end

%-----------------------------------------------
% convert units
%
% In mahowalds data, the unit of dust flux is 
%    [kg dust / m^2 / sec]
%
% the box model needs iron flux in 
%    [micromol Fe / m^2 / yr]
%
% (iron unit in model is [nmol/L] = [micromol/m^3])
%
% conversion:
% 1 kg dust = 1000 g dust
%           = 1000 * 0.035 g Fe 
%           = 1000 * 0.035 / 56 mol Fe
%           = 1000 * 0.035 / 56 * 1.0e6 micromol Fe
% in addition, 1 yr is 360*86400 seconds
%
% Note: We do not take into account the solubility of dust iron here;
% this is done internally in the model. However, a fixed percentage of
% iron in dust (3.5 percent) has been used here.
%-----------------------------------------------
unitconv = 3.5e7 / 56.0 * 360 * 86400; % to convert from kg dust/m^2/sec to 
                                       % micromol Fe/m^2/yr
dustflux = unitconv * mean(dust_clim_new,3);

%-----------------------------------------------
% calculate area of grid boxes for dust field
%-----------------------------------------------
dlon = mean(diff(lon)); % diff(lon) contains identical values
dlat0 = diff(lat(1:2));
latb = (lat(2:end) + lat(1:end-1))/2;
latb = [latb(1)-dlat0; latb; latb(end)+dlat0];
dlat = diff(latb);
earth_radius = 6371.0e3; 
nx = length(lon);
ny = length(lat);
rA = zeros(nx,ny);
for j=1:ny,
  dx = 2 * pi * earth_radius * dlon / 360 * cosd(lat(j));
  dy = 2 * pi * earth_radius * dlat(j) / 360;
  rA(:,j) = dx * dy;
end

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

box1_mask = (Indopac_mask_surf .* lat_mask_north);
box2_mask = (Indopac_mask_surf .* lat_mask_equat);
box3_mask = ((Indopac_mask_surf + Atlant_mask_surf) .* lat_mask_south);
box4_mask = (Atlant_mask_surf .* lat_mask_equat);
box5_mask = (Atlant_mask_surf .* lat_mask_north);

%------------------
% for each grid point of the dust field, find to which of the five 
% boxmodel-boxes it belongs; we do this by nearest neighbor interpolation
%------------------
boxnumbers = box1_mask + 2*box2_mask + 3*box3_mask + ...
    4*box4_mask + 5*box5_mask;
boxnumbers_dust = interp2(lat2,lon2',boxnumbers,lat,lon','nearest');

%------------------
% now integrate the dust flux over each of the five boxes
% and write out the result
%------------------

totdust = zeros(5,1);
for k = 1:5,
  ii = find(boxnumbers_dust==k);
  totdust(k) = sum( rA(ii) .* dustflux(ii) );
  fprintf('flux into box %i: %12.6e micromol/yr \n',k,totdust(k));
end
ii = find(boxnumbers_dust==0);
dust_on_land = sum( rA(ii) .* dustflux(ii) );
fprintf('total flux onto land:  %12.6e micromol/yr \n',dust_on_land);
fprintf('total flux into ocean: %12.6e micromol/yr \n',sum(totdust));


