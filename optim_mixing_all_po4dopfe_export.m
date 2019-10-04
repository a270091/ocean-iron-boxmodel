% optimizes export production to fit the PO4 distribution in the boxes

pvec_ini = zeros(7,1);
% mixing rates (in Sv) between several surface boxes
pvec_ini(1) = 40;
pvec_ini(2) = 40;
pvec_ini(3) = 17;
pvec_ini(4) = 10;
pvec_ini(5) = 5;
pvec_ini(6) = 20;
pvec_ini(7) = 10;

options = optimset('Display','iter','TolX',0.1,'TolFun',0.02);
[f_ini,po4_init0,po4_end0] = costf_mix_boxmodel_po4dopfe_export(pvec_ini);
pvec = fminsearch(@costf_mix_boxmodel_po4dopfe_export, pvec_ini,options);
rplot = 1;
[f_end,po4_init,po4_end] = costf_mix_boxmodel_po4(pvec,rplot);

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

