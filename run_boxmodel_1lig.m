% global fields for analysis
global params
global rhs feprime export

%-----------------------------
% initialize model parameters
%-----------------------------

modelrun = '1H';         % 1L or 1H or generic
if strcmp(modelrun,'1L')
    % settings corresponding to run 1L in Frontiers paper
    fprintf('Parameter setting: 1L\n');
    params = load('results/parameters_1l.mat');
elseif strcmp(modelrun,'1H')
    % settings corresponding to run 1L in Frontiers paper
    fprintf('Parameter setting: 1H\n');
    params = load('results/parameters_3l2.mat');
else
    % generic setting of parameters
    fprintf('Generic parameter setting\n');
    boxmodel_init_params()
end

% initial PO4 distribution
po4_init = params.po4init;
dop_init = zeros(size(po4_init));
fe_init = 0.6 + zeros(size(po4_init));
lig_init = 1.0 + zeros(size(po4_init));
conc_init = [po4_init;dop_init;fe_init;lig_init];

% integrate
tspan = (0:50:20000);
conc = ode23s(@boxmodel_dgl_po4dopfelig_export, tspan, conc_init);

% some diagnostics
totpo4 = params.volume' * conc.y(1:12,:);
totfe = params.volume' * conc.y(25:36,:);
c_export = export * params.redfield_c2p * 12 * 1.0e-18;

for k=1:12
  fprintf('Box %i: %s\n',k,params.long_names{k});
  fprintf('  DFe: %5.2f, Fe: %5.3f\n',conc.y(k+24,end), feprime(k));
  fprintf('  PO4: %5.2f, obs. PO4: %5.2f\n',conc.y(k,end),params.po4init(k));
  fprintf('  DOP: %5.2f\n',conc.y(k+12,end));
  fprintf('  Lig: %5.2f\n',conc.y(k+36,end));
  if k<=5,
    fprintf('  Export: %5.2f PgC/yr \n', c_export(k));
  end
  fprintf('\n');
end

% calculate some global numbers and print them

dfe_av = sum( conc.y(25:36,end) .* params.volume ) / sum(params.volume);
dfe_tot = sum( conc.y(25:36,end) .* params.volume ) * 1.0e-6; 
dust   = params.dust;
sol    = params.dustsol;
dust_fe_sol = sum( sol*dust ) * 1.0e-6;
sedfe  = params.sediment_fe;
sedfac = params.sed_fac;
sed_fe_sol = sum(sedfac * sedfe .* params.volume ) * 1.0e-6;
hydro  = params.hydrothermal;
hydfac = params.hydro_fac;
hyd_fe_sol = sum(hydfac * hydro .* params.volume ) * 1.0e-6;
residence_time = dfe_tot / (dust_fe_sol + sed_fe_sol + hyd_fe_sol);
fprintf('Total export: %5.2f PgC/yr \n', sum(c_export));
fprintf('Average dFe concentration: %5.2f nmol/L \n',dfe_av);
fprintf('Dust iron input: %8.2e mol/yr \n',dust_fe_sol);
fprintf('Sediment iron input: %8.2e mol/yr \n',sed_fe_sol);
fprintf('Hydrothermal iron input: %8.2e mol/yr \n',hyd_fe_sol);
fprintf('Residence time: %6.2f yr \n',residence_time);

% calculate some measures of goodness of model fit to data
ifit=1;
if (ifit==1)
    % model state
    finalstate = conc.y(:,end);
    po4mod = finalstate(1:12);
    femod  = finalstate(25:36);
    % observational data
    po4obs = params.po4init;
    read_geotraces_idp2;
    feobs  = femedian;
    % weight
    weight = ones(size(po4obs));
    [R,E,B,sf,sr] = calculate_taylor(po4obs,po4mod,weight);
    corr_po4 = R;
    rmse_po4 = E;
    bias_po4 = B;
    fprintf('\nModel assessment, PO4: \n');
    fprintf('  Corr: %f\n',corr_po4);
    fprintf('  RMSE: %f\n',rmse_po4);
    fprintf('  Bias: %f\n',bias_po4);
    % calculate for iron data
    [R,E,B,sf,sr] = calculate_taylor(feobs,femod,weight);
    corr_dfe = R;
    rmse_dfe = E;
    bias_dfe = B;
    mae_dfe  = mean( abs( feobs - femod ) );
    fprintf('\nModel assessment, Fe: \n');
    fprintf('  Corr: %f\n',corr_dfe);
    fprintf('  RMSE: %f\n',rmse_dfe);
    fprintf('  Bias: %f\n',bias_dfe);
    fprintf('  MAE:  %f\n',mae_dfe);
end

ifetable=0;
if (ifetable==1)
    % write out a table of iron sources and sinks for the final steady state
    finalstate = conc.y(:,end);
    dcdt = boxmodel_dgl_po4dopfelig_export(0,finalstate);
    % calculate fluxes (in 10^6 mol/yr) from rhs (calculated in the DGL field)
    dust  = rhs.dust  .* params.volume * 1.0e-12;
    sedfe = rhs.sedfe .* params.volume * 1.0e-12;
    hydro = rhs.hydro .* params.volume * 1.0e-12;
    bio   = (rhs.uptake + rhs.remin) .* params.volume * 1.0e-12;
    scav  = rhs.scav .* params.volume * 1.0e-12;
    fprintf('\niron fluxes: \n')
    for k=1:12
        fprintf('%7.1f %7.1f %7.1f %7.1f %7.1f\n',dust(k),sedfe(k),...
                hydro(k),bio(k),scav(k)) 
    end
end

iligtable=0;
if (iligtable==1)
    % write out a table of iron sources and sinks for the final steady state
    finalstate = conc.y(:,end);
    dcdt = boxmodel_dgl_po4dopfelig_export(0,finalstate);
    % calculate fluxes (in 10^6 mol/yr) from rhs (calculated in the DGL field)
    ligs1 = rhs.lig1  .* params.volume * 1.0e-12;
    ligs2 = rhs.lig2 .* params.volume * 1.0e-12;
    ligs3 = rhs.lig3 .* params.volume * 1.0e-12;
    fprintf('\nligand fluxes: \n')
    for k=1:12
        fprintf('%7.1f %7.1f %7.1f\n',ligs1(k),ligs2(k),...
                ligs3(k)) 
    end
end

iprint=0;
if (iprint==1)
    % save final concentrations as a table
    % (for later plotting of equilibrium ligands)
    finalstate = conc.y(:,end);
    fid = fopen('equil_po4dopfelig_export.dat','w');
    for k=1:12
        fprintf(fid,'%8.4f %8.3f %8.4f %8.4f\n', finalstate(k),...
                finalstate(k+12),finalstate(k+24),finalstate(k+36));
    end
    fclose(fid);
end

do_plot=0;

if (do_plot),
% plot phosphorus vs data
figure(1)
po4_end = conc.y(1:12,end);
hh = plot(po4_init, po4_end,'x');
hold on 
set(hh,'LineWidth',2,'MarkerSize',12,'Color',[0 0.6 0.3]);
set(gca,'PlotBoxAspectRatio',[1 1 1],'FontSize',12);
xlabel('WOA-derived PO_4');
ylabel('model PO_4');
yl = get(gca,'YLim');
hd = plot(yl,yl,'k');
print('PO4_boxmodel_po4dopfelig_sed.png','-dpng')

% plot phosphorus for each box, with data
%figure
%h = plot(conc.y(1:12,end),'kx');
%set(h,'MarkerSize',10, 'LineWidth',2);
%hold on
% hd = plot(params.po4init,'rx');
% set(hd,'MarkerSize',10, 'LineWidth',2);
% ylabel('PO_4 [\mumol L^{-1}]');
% set(gca,'XTick',(1:12),'XTickLabel',params.names,'XTickLabelRotation',45.0);
% set(gca,'FontSize',12,'XLim',[0.5,12.5],'YLim',[0 3.5]);

% read in iron data and make a plot of data with observations and model
figure(2)
plot_geotraces_idp2;
hold on
h = plot(conc.y(25:36,end),'kx');
set(h,'MarkerSize',10, 'LineWidth',2);
ylabel('dFe [nmol L^{-1}]');
set(gca,'XTick',(1:12),'XTickLabel',params.names,'XTickLabelRotation',45.0);
set(gca,'FontSize',12,'XLim',[0.5,12.5],'YLim',[0 2]);
print('fe_vs_data_po4dopdfelig_sed.png','-dpng');

end
