%---------------------------------------------------------------------------------
% read output of feedback analysis for a combined plot for the different model versions
%---------------------------------------------------------------------------------

results1L = load('results/feedback_1L_dust.mat');
results1H = load('results/feedback_1H_dust.mat');
results2L = load('results/feedback_2L_dust.mat');
results2H = load('results/feedback_2H_dust.mat');

%---------------------------------------------------------------------------------
% define colors for the 'Bright' color scheme (https://personal.sron.nl/~pault/)
%---------------------------------------------------------------------------------
bpurple = [170, 51,119]/256;
bred    = [238,102,119]/256;
byellow = [204,187, 68]/256;
bgreen  = [ 34,136, 51]/256;
bcyan   = [102,204,238]/256;
bblue   = [ 86,119,170]/256;

%---------------------------------------------------------------------------------
% calculate feedback factors for runs 1L and 1H
%---------------------------------------------------------------------------------

ik0 = find(results1L.percentage==1);
f_export = (results1L.tot_export_run2(ik0+1) - results1L.tot_export_run2(ik0-1)) / ...
    (results1L.tot_export_run1(ik0+1) - results1L.tot_export_run1(ik0-1));
f_surf_fe = (results1L.surf_fe_run2(ik0+1) - results1L.surf_fe_run2(ik0-1)) / ...
    (results1L.surf_fe_run1(ik0+1) - results1L.surf_fe_run1(ik0-1));
f_av_fe = (results1L.av_fe_run2(ik0+1) - results1L.av_fe_run2(ik0-1)) / ...
    (results1L.av_fe_run1(ik0+1) - results1L.av_fe_run1(ik0-1));
f_SO_fe = (results1L.SO_fe_run2(ik0+1) - results1L.SO_fe_run2(ik0-1)) / ...
    (results1L.SO_fe_run1(ik0+1) - results1L.SO_fe_run1(ik0-1));

fprintf('feedback factors and gains 1L \n')
fprintf('dFe (average): f =%5.2f g =%5.2f\n',f_av_fe,((f_av_fe-1)/f_av_fe))
fprintf('dFe (surface): f =%5.2f g =%5.2f\n',f_surf_fe,((f_surf_fe-1)/f_surf_fe))
fprintf('dFe (S.Ocean): f =%5.2f g =%5.2f\n',f_SO_fe,((f_SO_fe-1)/f_SO_fe))
fprintf('Export       : f =%5.2f g =%5.2f\n',f_export,((f_export-1)/f_export))

ik0 = find(results1H.percentage==1);
f_export = (results1H.tot_export_run2(ik0+1) - results1H.tot_export_run2(ik0-1)) / ...
    (results1H.tot_export_run1(ik0+1) - results1H.tot_export_run1(ik0-1));
f_surf_fe = (results1H.surf_fe_run2(ik0+1) - results1H.surf_fe_run2(ik0-1)) / ...
    (results1H.surf_fe_run1(ik0+1) - results1H.surf_fe_run1(ik0-1));
f_av_fe = (results1H.av_fe_run2(ik0+1) - results1H.av_fe_run2(ik0-1)) / ...
    (results1H.av_fe_run1(ik0+1) - results1H.av_fe_run1(ik0-1));
f_SO_fe = (results1H.SO_fe_run2(ik0+1) - results1H.SO_fe_run2(ik0-1)) / ...
    (results1H.SO_fe_run1(ik0+1) - results1H.SO_fe_run1(ik0-1));

fprintf('feedback factors and gains 1H\n')
fprintf('dFe (average): f =%5.2f g =%5.2f\n',f_av_fe,((f_av_fe-1)/f_av_fe))
fprintf('dFe (surface): f =%5.2f g =%5.2f\n',f_surf_fe,((f_surf_fe-1)/f_surf_fe))
fprintf('dFe (S.Ocean): f =%5.2f g =%5.2f\n',f_SO_fe,((f_SO_fe-1)/f_SO_fe))
fprintf('Export       : f =%5.2f g =%5.2f\n',f_export,((f_export-1)/f_export))

%---------------------------------------------------------------------------------
% calculate feedback factors for runs 2L and 2H
%---------------------------------------------------------------------------------

ik0 = find(results2L.percentage==1);
f_export = (results2L.tot_export_run2(ik0+1) - results2L.tot_export_run2(ik0-1)) / ...
    (results2L.tot_export_run1(ik0+1) - results2L.tot_export_run1(ik0-1));
f_surf_fe = (results2L.surf_fe_run2(ik0+1) - results2L.surf_fe_run2(ik0-1)) / ...
    (results2L.surf_fe_run1(ik0+1) - results2L.surf_fe_run1(ik0-1));
f_av_fe = (results2L.av_fe_run2(ik0+1) - results2L.av_fe_run2(ik0-1)) / ...
    (results2L.av_fe_run1(ik0+1) - results2L.av_fe_run1(ik0-1));
f_SO_fe = (results2L.SO_fe_run2(ik0+1) - results2L.SO_fe_run2(ik0-1)) / ...
    (results2L.SO_fe_run1(ik0+1) - results2L.SO_fe_run1(ik0-1));

fprintf('feedback factors and gains 2L \n')
fprintf('dFe (average): f =%5.2f g =%5.2f\n',f_av_fe,((f_av_fe-1)/f_av_fe))
fprintf('dFe (surface): f =%5.2f g =%5.2f\n',f_surf_fe,((f_surf_fe-1)/f_surf_fe))
fprintf('dFe (S.Ocean): f =%5.2f g =%5.2f\n',f_SO_fe,((f_SO_fe-1)/f_SO_fe))
fprintf('Export       : f =%5.2f g =%5.2f\n',f_export,((f_export-1)/f_export))

ik0 = find(results2H.percentage==1);
f_export = (results2H.tot_export_run2(ik0+1) - results2H.tot_export_run2(ik0-1)) / ...
    (results2H.tot_export_run1(ik0+1) - results2H.tot_export_run1(ik0-1));
f_surf_fe = (results2H.surf_fe_run2(ik0+1) - results2H.surf_fe_run2(ik0-1)) / ...
    (results2H.surf_fe_run1(ik0+1) - results2H.surf_fe_run1(ik0-1));
f_av_fe = (results2H.av_fe_run2(ik0+1) - results2H.av_fe_run2(ik0-1)) / ...
    (results2H.av_fe_run1(ik0+1) - results2H.av_fe_run1(ik0-1));
f_SO_fe = (results2H.SO_fe_run2(ik0+1) - results2H.SO_fe_run2(ik0-1)) / ...
    (results2H.SO_fe_run1(ik0+1) - results2H.SO_fe_run1(ik0-1));

fprintf('feedback factors and gains 2H\n')
fprintf('dFe (average): f =%5.2f g =%5.2f\n',f_av_fe,((f_av_fe-1)/f_av_fe))
fprintf('dFe (surface): f =%5.2f g =%5.2f\n',f_surf_fe,((f_surf_fe-1)/f_surf_fe))
fprintf('dFe (S.Ocean): f =%5.2f g =%5.2f\n',f_SO_fe,((f_SO_fe-1)/f_SO_fe))
fprintf('Export       : f =%5.2f g =%5.2f\n',f_export,((f_export-1)/f_export))

%---------------------------------------------------------------------------------
% make a plot of the linearized and the nonlinear reaction for the cases with one ligand
%---------------------------------------------------------------------------------

figure(1)
h1 = plot(results1L.percentage*100, results1L.tot_export_run1);
set(h1,'LineWidth',2,'Color',bblue,'LineStyle','--')
hold on
h2 = plot(results1L.percentage*100, results1L.tot_export_run2);
set(h2,'LineWidth',2,'Color',bblue)
h3 = plot(results1H.percentage*100, results1H.tot_export_run1);
set(h3,'LineWidth',2,'Color',bgreen,'LineStyle','--')
h4 = plot(results1H.percentage*100, results1H.tot_export_run2);
set(h4,'LineWidth',2,'Color',bgreen)
set(gca,'FontSize',12)
xlabel('change in total dust (%)')
ylabel('total export production (PgC yr^{-1})')
xl  = get(gca,'XLim');
yl  = get(gca,'YLim');
xpos = xl(1) + 0.85*(xl(2) - xl(1));
ypos = yl(1) + 0.1 *(yl(2) - yl(1));
ht = text(xpos,ypos,'d');
set(ht,'Fontweight','b','FontSize',30)
print('feedback_1lig_dust_export.png','-dpng')

figure(2)
h1 = plot(results1L.percentage*100, results1L.surf_fe_run1);
set(h1,'LineWidth',2,'Color',bblue,'LineStyle','--')
hold on
h2 = plot(results1L.percentage*100, results1L.surf_fe_run2);
set(h2,'LineWidth',2,'Color',bblue)
h3 = plot(results1H.percentage*100, results1H.surf_fe_run1);
set(h3,'LineWidth',2,'Color',bgreen,'LineStyle','--')
h4 = plot(results1H.percentage*100, results1H.surf_fe_run2);
set(h4,'LineWidth',2,'Color',bgreen)
set(gca,'FontSize',12)
xlabel('change in total dust (%)')
ylabel('global average surface dFe (nmol L^{-1})')
xl  = get(gca,'XLim');
yl  = get(gca,'YLim');
xpos = xl(1) + 0.85*(xl(2) - xl(1));
ypos = yl(1) + 0.1 *(yl(2) - yl(1));
ht = text(xpos,ypos,'b');
set(ht,'Fontweight','b','FontSize',30)
print('feedback_1lig_dust_surf_fe.png','-dpng')

figure(3)
h1 = plot(results1L.percentage*100, results1L.av_fe_run1);
set(h1,'LineWidth',2,'Color',bblue,'LineStyle','--')
hold on
h2 = plot(results1L.percentage*100, results1L.av_fe_run2);
set(h2,'LineWidth',2,'Color',bblue)
h3 = plot(results1H.percentage*100, results1H.av_fe_run1);
set(h3,'LineWidth',2,'Color',bgreen,'LineStyle','--')
h4 = plot(results1H.percentage*100, results1H.av_fe_run2);
set(h4,'LineWidth',2,'Color',bgreen)
set(gca,'FontSize',12)
xlabel('change in total dust (%)')
ylabel('global average dFe (nmol L^{-1})')
xl  = get(gca,'XLim');
yl  = get(gca,'YLim');
xpos = xl(1) + 0.85*(xl(2) - xl(1));
ypos = yl(1) + 0.1 *(yl(2) - yl(1));
ht = text(xpos,ypos,'a');
set(ht,'Fontweight','b','FontSize',30)
print('feedback_1lig_dust_ave_fe.png','-dpng')

figure(4)
h1 = plot(results1L.percentage*100, results1L.SO_fe_run1);
set(h1,'LineWidth',2,'Color',bblue,'LineStyle','--')
hold on
h2 = plot(results1L.percentage*100, results1L.SO_fe_run2);
set(h2,'LineWidth',2,'Color',bblue)
h3 = plot(results1H.percentage*100, results1H.SO_fe_run1);
set(h3,'LineWidth',2,'Color',bgreen,'LineStyle','--')
h4 = plot(results1H.percentage*100, results1H.SO_fe_run2);
set(h4,'LineWidth',2,'Color',bgreen)
set(gca,'FontSize',12)
xlabel('change in total dust (%)')
ylabel('Southern Ocean dFe (nmol L^{-1})')
xl  = get(gca,'XLim');
yl  = get(gca,'YLim');
xpos = xl(1) + 0.85*(xl(2) - xl(1));
ypos = yl(1) + 0.1 *(yl(2) - yl(1));
ht = text(xpos,ypos,'c');
set(ht,'Fontweight','b','FontSize',30)
print('feedback_1lig_dust_SO_fe.png','-dpng')

%---------------------------------------------------------------------------------
% now the linearized and the nonlinear reaction for the cases with two ligands
%---------------------------------------------------------------------------------

figure(5)
h1 = plot(results2L.percentage*100, results2L.tot_export_run1);
set(h1,'LineWidth',2,'Color',bcyan,'LineStyle','--')
hold on
h2 = plot(results2L.percentage*100, results2L.tot_export_run2);
set(h2,'LineWidth',2,'Color',bcyan)
h3 = plot(results2H.percentage*100, results2H.tot_export_run1);
set(h3,'LineWidth',2,'Color',byellow,'LineStyle','--')
h4 = plot(results2H.percentage*100, results2H.tot_export_run2);
set(h4,'LineWidth',2,'Color',byellow)
set(gca,'FontSize',12)
xlabel('change in total dust (%)')
ylabel('total export production (PgC yr^{-1})')
xl  = get(gca,'XLim');
yl  = get(gca,'YLim');
xpos = xl(1) + 0.85*(xl(2) - xl(1));
ypos = yl(1) + 0.1 *(yl(2) - yl(1));
ht = text(xpos,ypos,'d');
set(ht,'Fontweight','b','FontSize',30)
print('feedback_2lig_dust_export.png','-dpng')

figure(6)
h1 = plot(results2L.percentage*100, results2L.surf_fe_run1);
set(h1,'LineWidth',2,'Color',bcyan,'LineStyle','--')
hold on
h2 = plot(results2L.percentage*100, results2L.surf_fe_run2);
set(h2,'LineWidth',2,'Color',bcyan)
h3 = plot(results2H.percentage*100, results2H.surf_fe_run1);
set(h3,'LineWidth',2,'Color',byellow,'LineStyle','--')
h4 = plot(results2H.percentage*100, results2H.surf_fe_run2);
set(h4,'LineWidth',2,'Color',byellow)
set(gca,'FontSize',12)
xlabel('change in total dust (%)')
ylabel('global average surface dFe (nmol L^{-1})')
xl  = get(gca,'XLim');
yl  = get(gca,'YLim');
xpos = xl(1) + 0.85*(xl(2) - xl(1));
ypos = yl(1) + 0.1 *(yl(2) - yl(1));
ht = text(xpos,ypos,'b');
set(ht,'Fontweight','b','FontSize',30)
print('feedback_2lig_dust_surf_fe.png','-dpng')

figure(7)
h1 = plot(results2L.percentage*100, results2L.av_fe_run1);
set(h1,'LineWidth',2,'Color',bcyan,'LineStyle','--')
hold on
h2 = plot(results2L.percentage*100, results2L.av_fe_run2);
set(h2,'LineWidth',2,'Color',bcyan)
h3 = plot(results2H.percentage*100, results2H.av_fe_run1);
set(h3,'LineWidth',2,'Color',byellow,'LineStyle','--')
h4 = plot(results2H.percentage*100, results2H.av_fe_run2);
set(h4,'LineWidth',2,'Color',byellow)
set(gca,'FontSize',12)
xlabel('change in total dust (%)')
ylabel('global average dFe (nmol L^{-1})')
xl  = get(gca,'XLim');
yl  = get(gca,'YLim');
xpos = xl(1) + 0.85*(xl(2) - xl(1));
ypos = yl(1) + 0.1 *(yl(2) - yl(1));
ht = text(xpos,ypos,'a');
set(ht,'Fontweight','b','FontSize',30)
print('feedback_2lig_dust_ave_fe.png','-dpng')

figure(8)
h1 = plot(results2L.percentage*100, results2L.SO_fe_run1);
set(h1,'LineWidth',2,'Color',bcyan,'LineStyle','--')
hold on
h2 = plot(results2L.percentage*100, results2L.SO_fe_run2);
set(h2,'LineWidth',2,'Color',bcyan)
h3 = plot(results2H.percentage*100, results2H.SO_fe_run1);
set(h3,'LineWidth',2,'Color',byellow,'LineStyle','--')
h4 = plot(results2H.percentage*100, results2H.SO_fe_run2);
set(h4,'LineWidth',2,'Color',byellow)
set(gca,'FontSize',12)
xlabel('change in total dust (%)')
ylabel('Southern Ocean dFe (nmol L^{-1})')
xl  = get(gca,'XLim');
yl  = get(gca,'YLim');
xpos = xl(1) + 0.85*(xl(2) - xl(1));
ypos = yl(1) + 0.1 *(yl(2) - yl(1));
ht = text(xpos,ypos,'c');
set(ht,'Fontweight','b','FontSize',30)
print('feedback_2lig_dust_SO_fe.png','-dpng')

