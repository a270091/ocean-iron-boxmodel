function S = read_one_netcdf(fname)
% testscript to read one chlorophyll file

nc = netcdf.open(fname,'NC_NOWRITE');

% get basic information about netcdf file
[ndims nvars natts recdim] = netcdf.inq(nc);

varlist = {};
for k=0:nvars-1,
  varlist{k+1} = netcdf.inqVar(nc,k);
end

S.dimensions = [];
dimlist = {};
for k=0:ndims-1,
  [dimname,dimsize] = netcdf.inqDim(nc,k);
  S.dimensions.(dimname) = dimsize;
end
 
% Attributes for structure
S.attributes.global=ncgetatt(nc,'global');
      
% Read variable data
for ivar=1:size(varlist,2)
    
  cvar=char(varlist{ivar});
  varid=ncfindvarid(nc,cvar);
  if isempty(varid)
    disp(['No such variable ''',cvar,''' in MNC file ',fname]);
    continue
  end
    
  [varname,xtype,dimids,natts] = netcdf.inqVar(nc,varid);
  vf = double(netcdf.getVar(nc,varid));

  S.(cvar)=vf;
  S.attributes.(cvar)=ncgetatt(nc,cvar);
%  S.dimensions.(cvar)= ncgetdims(nc,cvar);
  S.var_dims.(cvar)= ncgetdims(nc,cvar);
%  S.dimensions.(cvar) = 
 
end

netcdf.close(nc);
if isempty(S)
  error('Something didn''t work!!!');
end

return
  
function A = ncgetatt(nc,varname)
% get all attributes and put them into a struct
  
% 1. get global properties of file
[ndims nvars natts dimm] = netcdf.inq(nc);

% get variable ID and properties
if strcmp('global',varname)
  % "variable" is global
  varid = netcdf.getConstant('NC_GLOBAL');
else
  % find variable ID and properties
  varid = ncfindvarid(nc,varname);
  if ~isempty(varid)
    [varn,xtype,dimids,natts] = netcdf.inqVar(nc,varid);
  else
    warning(sprintf('variable %s not found',varname))
  end
end

if natts >= 1
  for k=0:natts-1
    attn = netcdf.inqAttName(nc,varid,k);
    [xtype attlen]=netcdf.inqAtt(nc,varid,attn);
    attval = netcdf.getAtt(nc,varid,attn);
    if ~ischar(attval)
      attval = double(attval);
    end
    if strcmp(attn(1),'_'),
      attn = attn(2:end);
    end
    A.(char(attn))=attval;
  end
else
  A = 'none';
end
  
return
  
function [dimvec,sizes] = ncgetdims(nc,varname)
% get all dimensions that belong to a variable as a cell array
% of strings and get the size of the variable (the latter is
% not really new information)
  
varid = ncfindvarid(nc,varname);
if ~isempty(varid)
  [varn,xtype,dimids,natts] = netcdf.inqVar(nc,varid);
else
  warning(sprintf('variable %s not found',varname))
end
  
nvardim = length(dimids);
dimvec = {};
sizes = [];
  
if nvardim > 0
  for k=1:nvardim,
    [dimname,dimlength] = netcdf.inqDim(nc,dimids(k));
    dimvec(k) = {dimname};
    sizes(k) = dimlength;
  end
end

return

function varid = ncfindvarid(nc,varname)

[ndims nvars natts dimm] = netcdf.inq(nc);
varid=[];
for k=0:nvars-1
  if strcmp(netcdf.inqVar(nc,k),varname);
    varid = netcdf.inqVarID(nc,varname);
    break
  end
end

return
