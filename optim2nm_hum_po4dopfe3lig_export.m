% optimizes parameters for the ligand dynamics to obtain a better
% fit with iron data

%----
% first step: read iron data and create median values for the boxes
%----
sort_fe_data_into_boxes

global fe_data pvec_dimensional
fe_data = femedian

%----
% initial guess for parameters (for normalization, we define the optimized
% parameters as simple factors multiplying the original dimensional parameters)
%----

pvec_ini = ones(4,1);
pvec_dimensional(1) = 0.3;      % params.hum
pvec_dimensional(2) = 2.5e-4 * 116; % Lig to P ratio in POC remineralization
                                    % (Lig:C ratio * Redfield C:P)
pvec_dimensional(3) = 1.8 * 2.5e-4 * 116 * 0.67;   % params.rlig2p2
pvec_dimensional(4) = 0.5e-3; % params.ligrem
pvec = pvec_ini;

[f_ini,dfe_ini] = costf2_hum_boxmodel_po4dopfe3lig_export(pvec_ini);

%----
% search for better parameters. Unlike the other optimization, 
% here we do a coarse latin hypecube followed by Nelder-Mead
%----
options = optimset('Display','iter','TolX',0.02,'TolFun',0.02);

ni = 3;
pvec_all = zeros(ni^3,4);
cost_all = zeros(ni^3,1);
dn = 2. / (ni+1);
nn = 0;
for k=1:ni
  pvec_ini(1) = dn*k;
  for l=1:ni,
    pvec_ini(2) = dn*l;
    for m=1:ni,
       pvec_ini(3) = dn*m;
       [pvec,cost] = fminsearch(@costf2_hum_boxmodel_po4dopfe3lig_export, pvec_ini,options);
       nn = nn+1;
       pvec_all(nn,:) = pvec(:);
       cost_all(nn) = cost;
    end
  end
end

%

% [f_end,dfe_end] = costf2_hum_boxmodel_po4dopfe3lig_export(pvec);
