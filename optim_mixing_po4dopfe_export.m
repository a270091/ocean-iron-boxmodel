% optimizes export production to fit the PO4 distribution in the boxes

pvec_ini = zeros(7,1);
% mixing rates (in Sv) between several surface boxes
pvec_ini(1) = 40;
pvec_ini(2) = 7;
pvec_ini(3) = 50;
pvec_ini(4) = 25;
pvec_ini(5) = 8;
pvec_ini(6) = 0;
pvec_ini(7) = 40;
pvec = pvec_ini;

% options = optimset('Display','iter','TolX',0.1,'TolFun',0.02);
[f_ini,po4_init0,po4_end0] = costf_mix_boxmodel_po4dopfe_export(pvec_ini);
%pvec = fminsearch(@costf_mix_boxmodel_po4, pvec_ini,options);

mixmin = 00;
mixmax = 70;
lmax = 41;
for l=1:lmax
  mixvec(l) = mixmin + (l-1)*(mixmax-mixmin)/(lmax-1);
  pvec(3) = mixmin + (l-1)*(mixmax-mixmin)/(lmax-1);
  c(l) = costf_mix_boxmodel_po4dopfe_export(pvec);
end
figure()
plot(mixvec,c)
[cmin,imin] = min(c);

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

