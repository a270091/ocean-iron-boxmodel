%------------------------------------------------------------
% Calculates the masks on the WOA grid for the boxes of a 
% 12-box model for the global biogeochemical cycles in the ocean.
%
% Then reads the GEOTRACES IDP data for Fe and sorts the iron 
% measurements into boxes to calculate median etc.  
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
% to calculate the volumes of specific water masses we need
% the areas and the vertical boundaries of the boxes 
%------------------------------------------------------------

r_earth = 6430.0e3; % mean earth radius

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

% create an 'index array', which contains for each WOA grid box
% to which box of the boxmodel it belongs
box_indices = 1*NPac_surf + 2*EqPac_surf + 3*SO_surf + ...
    4*EqAtl_surf + 5*NAtl_surf + 6*AAIW_Pac + 7*AAIW_Atl + ...
    8*Pac_deep + 9*SO_deep + 10*Atl_deep + ...
    11*Pac_bottom + 12*Atl_bottom;

%------------------------------------------------------------
% read data from GEOTRACES intermediate data product 2
%------------------------------------------------------------

% read data
location = '/Users/cvoelker/data/GEOTRACES/';
fname = 'GEOTRACES_IDP2017_v2_Discrete_Sample_Data_2c95e759_1.nc';
a = read_one_netcdf([location,fname]);

cruise = a.metavar1;
lon_fe      = a.longitude;
lat_fe      = a.latitude;
depth_fe    = a.var2;
ctdtemp  = a.var7;
ctdsalt  = a.var8;
fed      = a.var20;
fillval  = a.attributes.var20.FillValue;

% replace missing values
depth_fe(depth_fe==fillval) = NaN;
ctdtemp(ctdtemp==fillval)   = NaN;
ctdsalt(ctdsalt==fillval)   = NaN;
fed(fed==fillval)           = NaN;

nboxdata = zeros(12,1);
feboxdata = NaN*ones(12,1300);

nstations = length(lon_fe);
for l = 1:nstations
    ix = ceil(lon_fe(l));
    if ix==0, ix = 360; end
    jy = ceil(lat_fe(l))+90;
%     fprintf('lon %f lat %f %i %i\n', lon_fe(l), lat_fe(l),ix,jy)
    depthprof = depth_fe(:,l);
    feprof    = fed(:,l);
    % throw out missing values
    ii = find(~isnan(feprof) & ~isnan(depthprof));
    if isempty(ii), continue; end
    depthprof = depthprof(ii);
    feprof = feprof(ii);
%    fprintf('lon %f lat %f ndata %i\n', lon_fe(l), lat_fe(l),length(feprof))
    for k=1:length(feprof)
        % for each sample deoth find, n which WOA grid point it belongs
        sample_depth = depthprof(k);
        kk = sum(sample_depth>depth_bnds(1:33));
        whichbox = box_indices(ix,jy,kk);
	if (feprof(k)>10)
	  fprintf('Fe: %f, lat,lon,depth: %f %f %f, box %i\n',...
		  feprof(k),lon_fe(l), lat_fe(l), depthprof(k), whichbox)
	end
        if (whichbox>0)
            nboxdata(whichbox) = nboxdata(whichbox)+1;
            feboxdata(whichbox,nboxdata(whichbox)) = feprof(k);
        end
    end
end

% now calculate statistics (mean, median, quartiles, for each of the boxes

femedian=zeros(12,1);
fe1q = femedian;
fe3q = femedian;
for k=1:12
    ndata = nboxdata(k);
    fedata = feboxdata(k,1:ndata);
    femedian(k) = median(fedata);
    femin(k)    = min(fedata);
    femax(k)    = max(fedata);
    felower = fedata(fedata<=femedian(k));
    feupper = fedata(fedata>=femedian(k));
    fe1q(k) = median(felower);
    fe3q(k) = median(feupper);
    fprintf('Box %2i: N=%4i, 1stQ=%5.3f Median=%5.3f 3rdQ=%5.3f\n',...
        k,ndata,fe1q(k),femedian(k),fe3q(k));
end

% finally, save the output as an ASCII-file
fname = 'results/geotraces_idp2_dfe_data_boxed.dat';
fid = fopen(fname,'w');
for k=1:12
    fprintf(fid,'%2i %6i %8.4f %8.4f %8.4f %8.4f %8.4f\n',...
            k,nboxdata(k),femin(k),fe1q(k),femedian(k),fe3q(k),femax(k));
end
fclose(fid);

% and now do a plot
% myboxplot(femedian,fe1q,fe3q,femin,femax)



return



