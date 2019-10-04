function [f,po4_init,po4_final] = costf_mix_boxmodel_po4dopfe_export(pvec,rplot)
% f = costf_boxmodel_po4(pvec)

% initialize model parameters
global params
boxmodel_init_params()

if (nargin>0) 
% advection: first define a three-column matrix where the first 
% column gives the number of the box that the advection originates
% from, the second the box it goes into, and the third the flux
% in Sverdrup (10^6 m^3/s)

  adv_matrix = params.aux.adv_matrix;

% we also need some mixing between the surface boxes and the boxes below; 
% mixing is conceptualized as biderictional flux from box A to box B and
% vice versa.
  mix_matrix = params.aux.mix_matrix;
  mix_matrix(:,3) = pvec;

  tr = zeros(12,12);
  for k=1:21
    i = adv_matrix(k,1);
    j = adv_matrix(k,2);
    tr(i,i) = tr(i,i) - adv_matrix(k,3);
    tr(j,i) = tr(j,i) + adv_matrix(k,3);
  end

  mix = zeros(12,12);
  for k=1:7
    i = mix_matrix(k,1);
    j = mix_matrix(k,2);
    mix(i,i) = mix(i,i) - mix_matrix(k,3);
    mix(j,j) = mix(j,j) - mix_matrix(k,3);
    mix(j,i) = mix(j,i) + mix_matrix(k,3);
    mix(i,j) = mix(i,j) + mix_matrix(k,3);
  end
  % tr = tr - tr';

  params.advect = zeros(12,12);
  unitfac = 1.0e6 * 360 * 86400; % to convert volume transport from 10^6 m^3/s
                               % to m^3/yr
  for k=1:12
    params.advect(k,:) = unitfac * (tr(k,:) + mix(k,:))/ params.volume(k);
  end
end

% initial PO4 distribution
po4_init = params.po4init;
dop_init = zeros(size(po4_init));
fe_init = 0.6 + zeros(size(po4_init));
conc_init = [po4_init;dop_init;fe_init];

% integrate
tspan = (0:50:3000);
conc = ode23s(@boxmodel_dgl_po4dopfe_export, tspan, conc_init);

% calculate RMS difference between final and initial PO4
po4_final = conc.y(1:12,end);
diff = po4_init - po4_final;
% f = sqrt( sum( (diff.^2).*params.volume ) ./ sum(params.volume) );
f = sqrt( sum( (diff.^2) ) );
%f = sum( abs(diff ) );

% prevent taking too strong mixing
if (max(abs(pvec))>80) 
    f = f + max(abs(pvec))-80;
end
if (min(pvec)<0) 
    f = f + 1000;
end
if (min(po4_final)<0)
    f = f+1000;
end

if (nargin>1)
    if (rplot>0)
        plot(conc.x,conc.y);
    end
end

return
end
