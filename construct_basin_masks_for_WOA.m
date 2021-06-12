
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

return

