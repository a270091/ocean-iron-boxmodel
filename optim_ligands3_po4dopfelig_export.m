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
pvec_dimensional(1) = 2.5e-4 * 116;
pvec_dimensional(2) = 5.0e-4 * 116;
pvec_dimensional(3) = 0.5e-3;
pvec_dimensional(4) = 100;
pvec = pvec_ini;

[f_ini,dfe_ini] = costf2_ligands_boxmodel_po4dopfelig_export(pvec_ini);

%----
% search for better parameters. Unlike the other optimization, 
% here we just do a latin hypercube
%----
nk = 21;
misfit = zeros(nk,1);
p1 = zeros(nk,1);
for k=1:nk
  f = (k-1)/10 - 1.0;
  pvec(4) = 10.0^(f);
  p1(k) = pvec(4);
  misfit(k) = costf2_ligands_boxmodel_po4dopfelig_export(pvec);
end

plot(p1,misfit)

return

[f_end,dfe_end] = costf2_ligands_boxmodel_po4dopfelig_export(pvec);

return

% make plot of PO4 vs data for intial and best solution
pvec(3) = mixvec(imin);
[f_end,po4_init,po4_end] = costf_mix_boxmodel_po4dopfe_export(pvec);

figure()
hh = plot(po4_init, po4_end,'x');
set(hh,'LineWidth',2,'MarkerSize',12,'Color',[0 0.6 0.3]);
hold on
hh0 = plot(po4_init0, po4_end0,'x');
set(hh0,'LineWidth',2,'MarkerSize',12,'Color',[0.6 0.1 0.3]);
set(gca,'PlotBoxAspectRatio',[1 1 1],'FontSize',12);
xlabel('WOA-derived PO_4');
ylabel('model PO_4');
yl = get(gca,'YLim');
hd = plot(yl,yl,'k');

