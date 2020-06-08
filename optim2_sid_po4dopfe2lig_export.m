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
pvec_dimensional(1) = 6.0;      % params.beta
pvec_dimensional(2) = 0.1;      % params.KFe_bact
pvec_dimensional(3) = 1000.0;   % params.ksid
pvec_dimensional(4) = 5.0e-4 * 116;   % params.rlig2p2
pvec = pvec_ini;

[f_ini,dfe_ini] = costf_sid_boxmodel_po4dopfe2lig_export(pvec_ini);

%----
% search for better parameters. Unlike the other optimization, 
% here we just do a latin hypercube
%----

% here we change both the production rate and the ligand stability constant,
% assuming no production of weak ligand during phytoplankton growth
% (pvec(4)=0)
pvec(4) = 0.0;
nk = 11;
nm = 11;
misfit = zeros(nk,nm);
p1 = zeros(nk,1);
p1max = 20.0;
p2 = zeros(nk,1);
p2max = 0.25;
for k=1:nk
  pvec(1) = (k-1)*p1max/(nk-1);
  p1(k) = pvec(1);
  for m=1:nk
    fprintf('%i %i\n',k,m)
    pvec(2) = (m-1)/(nm-1)*p2max;
    p2(m) = pvec(4);
    misfit(k,m) = costf_sid_boxmodel_po4dopfe2lig_export(pvec);
  end
end

figure;
plot(p1,misfit(:,1));
hold on
plot(p1,misfit(:,6),'r');
plot(p1,misfit(:,11),'g');

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

