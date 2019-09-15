% optimizes export production to fit the PO4 distribution in the boxes

pvec_ini = zeros(5,1);
% export in PgC/yr
pvec_ini(1) = 0.99; % 1.87;
pvec_ini(2) = 1.88; % 2.83;
pvec_ini(3) = 2.38; % 2.66;
pvec_ini(4) = 0.72; % 1.20;
pvec_ini(5) = 1.45; % 1.53;

options = optimset('Display','iter','TolX',0.01,'TolFun',0.01);
[f_ini, po4_init, po4_end] = costf_exp_boxmodel_po4(pvec_ini);
pvec = fminsearch(@costf_exp_boxmodel_po4, pvec_ini,options);
rplot=1;
[f_end, po4_init0, po4_end0] = costf_exp_boxmodel_po4(pvec,rplot);

figure()
hh = plot(po4_init, po4_end,'x');
set(hh,'LineWidth',2,'MarkerSize',12,'Color',[0 0.6 0.3]);
hold on
hh0 = plot(po4_init0, po4_end0,'x');
set(hh0,'LineWidth',2,'MarkerSize',12,'Color',[0.6 0.1 0.3]);
set(gca,'PlotBoxAspectRatio',[1 1 1],'FontSize',12);
xlabel('WOA-derived PO$_4$');
ylabel('model PO$_4$');
yl = get(gca,'YLim');
hd = plot(yl,yl,'k');

