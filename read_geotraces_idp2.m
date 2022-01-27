%------------------------------------------------------------
% read data from GEOTRACES intermediate data product 2, already 
% sorted and collected into boxes
%------------------------------------------------------------

nboxdata = zeros(12,1);
femedian = zeros(12,1);
fe1q     = femedian;
fe3q     = femedian;
femin    = femedian;
femax    = femedian;

% finally, save the output as an ASCII-file
fname = 'results/geotraces_idp2_dfe_data_boxed.dat';
fid = fopen(fname,'r');
for k=1:12
    a = str2num(fgetl(fid));
    nboxdata(k) = a(2);
    femin(k)    = a(3);
    fe1q(k)     = a(4);
    femedian(k) = a(5);
    fe3q(k)     = a(6);
    femax(k)    = a(7);
end
fclose(fid);

return



