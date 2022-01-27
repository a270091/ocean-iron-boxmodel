% varies total Fe sources to the ocean and Fe scavenging rate over a
% pre-defined range, to study the sensitivity of iron concentrations
% to these factors

fes_min = 0.5;
fes_max = 10.0;
n_fes = 20;
sca_min = 1.0;
sca_max = 50.0;
n_sca = 50;

% initialize model parameters
global params
boxmodel_init_params()

% initial PO4 distribution
po4_init = params.po4init;
dop_init = zeros(size(po4_init));
fe_init = 0.6 + zeros(size(po4_init));
conc_init = [po4_init;dop_init;fe_init];

% save 'reference' values of Fe sources and scavenging rate
scav_ref  = params.kscav;
dust_ref  = params.dust;
hydro_ref = params.hydrothermal;
sed_ref   = params.sediment_fe;

for nx = 1:n_fes
    fes_fac(nx) = fes_min + (nx-1)*(fes_max - fes_min)/(n_fes-1);
    params.dust         = fes_fac(nx) * dust_ref;
    params.hydrothermal = fes_fac(nx) * hydro_ref;
    params.sediment_fe  = fes_fac(nx) * sed_ref;
    for ny = 1:n_sca
        sca_fac(ny) = sca_min + (nx-1)*(sca_max - sca_min)/(n_sca-1);
        params.kscav = sca_fac(ny) * scav_ref;
        % integrate
        tspan = (0:50:3000);
        conc = ode23s(@boxmodel_dgl_po4dopfe_export, tspan, conc_init);
        po4_final = conc.y(1:12,end);
        dop_final = conc.y(13:24,end);
        dfe_final = conc.y(25:36,end);
    end
end


