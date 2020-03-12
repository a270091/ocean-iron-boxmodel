% optimizes export production to fit the PO4 distribution in the boxes

pvec_ini(1) = 0.54;

% options = optimset('Display','iter','TolX',0.1,'TolFun',0.02);
[f_ini,po4_init0,po4_end0] = costf_remin_boxmodel_po4dopfe_export(pvec_ini);
%pvec = fminsearch(@costf_mix_boxmodel_po4, pvec_ini,options);

pvec = pvec_ini;
reminmin = 0.2;
reminmax = 0.6;
lmax = 31;
for l=1:lmax
  reminvec(l) = reminmin + (l-1)*(reminmax-reminmin)/(lmax-1);
  pvec(1) = reminvec(l);
  c(l) = costf_remin_boxmodel_po4dopfe_export(pvec);
end

% make a plot pf the costfunction
figure()
plot(reminvec,c)
[cmin,imin] = min(c);

% make plot of PO4 vs data for intial and best solution
pvec(1) = reminvec(imin);
[f_end,po4_init,po4_end] = costf_remin_boxmodel_po4dopfe_export(pvec);

fprintf('optimized remin: %5.3f\n', pvec)

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

