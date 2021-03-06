%---------------------------------------------------------------------------------
% Feedback analysis for a box model with variable ligands
% 
% The feedback analysis requires three different runs of the box model:
% - the model with prognostic ligands is run into equilibrium for standard 
%   external Fe input and other model parameters
% - the equilibrium ligand concentration from the first run is extracted and
%   a second model run is done where the external forcing is varied, 
%   keeping the ligand distribution fixed at the distribution from the first run
% - the model is run a third time, with varied external forcing, but this time 
%   allowing the ligands to change with repect to the first run. 
%
% We do this here in a loop, varying the percentage by which dust deposition is 
% changed from -50 percent to +50 percent, to be able to plot the response function
% (e.g. total surface iron) as a function of the percentage of perturbation.
%---------------------------------------------------------------------------------

% decide whether this program will produce plots or only output for some 
% other plotting routine
do_plot = 0; % 1 for plotting on

% global fields for analysis
global params
global rhs feprime export

% initialize model parameters
boxmodel_init_params()

% initial PO4, DOP, Fe and Ligand distribution
po4_init = params.po4init;
dop_init = zeros(size(po4_init));
fe_init = 0.6 + zeros(size(po4_init));
lig_init = 1.0 + zeros(size(po4_init));
sid_init = zeros(size(po4_init));
conc_init = [po4_init;dop_init;fe_init;lig_init;sid_init];

% Reset some parameters to better values for the siderophore.
% From the large number of different parameter sets, we use here two
% examples that are more or less equally good in terms of reproducing Fe
% observations. One has a siderophore lifetime of 200 years, the other of
% 20 years
% $$$ parchoice=2;
% $$$ if parchoice==1
% $$$   params.beta = 6.0 * 0.225;
% $$$   params.KFe_bact = 0.1 * 0.01;
% $$$   params.rlig2p2 = 2.5e-4 * 116 * 0.67;   % params.rlig2p2 = lig2p * dopfrac
% $$$   params.sidremin = 0.5*0.01;
% $$$ else
% $$$   params.beta = 6.0;
% $$$   params.KFe_bact = 0.1 * 0.01;
% $$$   params.rlig2p2 = 1.8 * 2.5e-4 * 116 * 0.67;   % params.rlig2p2 = lig2p * dopfrac
% $$$   params.sidremin = 0.5*0.1;
% $$$ end

% depending on which case we want to solve, overwrite the parameters
% with the one belonging to that case
model_run = '2H'; 
if (model_run=='2L'),
    params2 = load('results/parameters_2l2.mat');
    params = params2;
elseif (model_run=='2H'),
    params2 = load('results/parameters_3l1.mat');
    params = params2;
else
    fprintf('must be 1L or 1H \n')
end

%---------------------------------------------------------------------------------
% run0: first integration of the box model, with variable ligands
%---------------------------------------------------------------------------------

tspan = (0:50:5000);
conc = ode23s(@boxmodel_dgl_po4dopfe2lig_export, tspan, conc_init);

% final state of the model, and some diagnostics
conc_final_run0 = conc.y(:,end);
c_export_run0 = export * params.redfield_c2p * 12 * 1.0e-18;

%---------------------------------------------------------------------------------
% changing parameters. Here it may be interesting to check a bit which ones one 
% would change. In this case we change total dust deposition
%---------------------------------------------------------------------------------

dust_unperturbed = params.dust;

nk = 21;
percentage = zeros(nk,1);
tot_export_run1 = zeros(nk,1);
tot_export_run2 = zeros(nk,1);
surf_fe_run1    = zeros(nk,1);
surf_fe_run2    = zeros(nk,1);
av_fe_run1      = zeros(nk,1);
av_fe_run2      = zeros(nk,1);
SO_fe_run1      = zeros(nk,1);
SO_fe_run2      = zeros(nk,1);

for kl = 1:nk
  percentage(kl) = (kl-1)*0.1;
  fprintf('dust factor: %5.2f \n',percentage(kl))

  params.dust = dust_unperturbed*percentage(kl); % increase total dust by XX percent

%---------------------------------------------------------------------------------
% run1: box model, with ligands kept fixed at distributions from run0,
% i.e. without the feedback
%---------------------------------------------------------------------------------

  global ligfix sidfix
  ligfix = conc_final_run0(37:48);
  sidfix = conc_final_run0(49:60);

  conc_init_nofback = conc_init(1:36);
  tspan = (0:50:5000);
  conc = ode23s(@boxmodel_dgl_po4dopfe_2ligfix_export, tspan, conc_init_nofback);

  % final state of the model, and some diagnostics
  conc_final_run1 = conc.y(:,end);
  c_export_run1 = export * params.redfield_c2p * 12 * 1.0e-18;
  
  % define several quantities whose change we are interested in
  tot_export_run1(kl) = sum(c_export_run1);
  surf_fe_run1(kl) = sum(conc_final_run1(25:29) .* params.area) / sum(params.area);
  av_fe_run1(kl)   = sum(conc_final_run1(25:36) .* params.volume) / sum(params.volume);
  SO_fe_run1(kl)   = conc_final_run1(27);

%---------------------------------------------------------------------------------
% run2: box model, with ligands reacting to the changed parameters, i.e. 
% with the ligand feedback included
%---------------------------------------------------------------------------------

  tspan = (0:50:5000);
  conc = ode23s(@boxmodel_dgl_po4dopfe2lig_export, tspan, conc_init);

  % final state of the model, and some diagnostics
  conc_final_run2 = conc.y(:,end);
  c_export_run2 = export * params.redfield_c2p * 12 * 1.0e-18;

  % define several quantities whose change we are interested in
  tot_export_run2(kl) = sum(c_export_run2);
  surf_fe_run2(kl) = sum(conc_final_run2(25:29) .* params.area) / sum(params.area);
  av_fe_run2(kl)   = sum(conc_final_run2(25:36) .* params.volume) / sum(params.volume);
  SO_fe_run2(kl)   = conc_final_run2(27);

end

%---------------------------------------------------------------------------------
% save output for a possible combined plot for the different model versions
%---------------------------------------------------------------------------------
% put everything into a structure
results.percentage = percentage;
results.tot_export_run1 = tot_export_run1;
results.tot_export_run2 = tot_export_run2;
results.surf_fe_run1 = surf_fe_run1;
results.surf_fe_run2 = surf_fe_run2;
results.av_fe_run1 = av_fe_run1;
results.av_fe_run2 = av_fe_run2;
results.SO_fe_run1 = SO_fe_run1;
results.SO_fe_run2 = SO_fe_run2;

if (model_run=='2L'),
    save('feedback_2L_dust.mat','-struct','results')
elseif (model_run=='2H'),
    save('feedback_2H_dust.mat','-struct','results')
end

%---------------------------------------------------------------------------------
% calculate feedback factors
%---------------------------------------------------------------------------------

ik0 = find(percentage==1);
f_export = (tot_export_run2(ik0+1) - tot_export_run2(ik0-1)) / ...
    (tot_export_run1(ik0+1) - tot_export_run1(ik0-1));
f_surf_fe = (surf_fe_run2(ik0+1) - surf_fe_run2(ik0-1)) / ...
    (surf_fe_run1(ik0+1) - surf_fe_run1(ik0-1));
f_av_fe = (av_fe_run2(ik0+1) - av_fe_run2(ik0-1)) / ...
    (av_fe_run1(ik0+1) - av_fe_run1(ik0-1));
f_SO_fe = (SO_fe_run2(ik0+1) - SO_fe_run2(ik0-1)) / ...
    (SO_fe_run1(ik0+1) - SO_fe_run1(ik0-1));

fprintf('feedback factors and gains\n')
fprintf('dFe (average): f =%5.2f g =%5.2f\n',f_av_fe,((f_av_fe-1)/f_av_fe))
fprintf('dFe (surface): f =%5.2f g =%5.2f\n',f_surf_fe,((f_surf_fe-1)/f_surf_fe))
fprintf('dFe (S.Ocean): f =%5.2f g =%5.2f\n',f_SO_fe,((f_SO_fe-1)/f_SO_fe))
fprintf('Export       : f =%5.2f g =%5.2f\n',f_export,((f_export-1)/f_export))

if ~do_plot, return; end

%---------------------------------------------------------------------------------
% make a plot of the linearized and the nonlinear reaction
%---------------------------------------------------------------------------------

figure(1)
h1 = plot(percentage*100, tot_export_run1);
set(h1,'LineWidth',2)
hold on
h2 = plot(percentage*100, tot_export_run2);
set(h2,'LineWidth',2)
set(gca,'FontSize',12)
xlabel('change in total dust (%)')
ylabel('total export production (PgC yr^{-1})')
xl  = get(gca,'XLim');
yl  = get(gca,'YLim');
xpos = xl(1) + 0.85*(xl(2) - xl(1));
ypos = yl(1) + 0.1 *(yl(2) - yl(1));
ht = text(xpos,ypos,'d');
set(ht,'Fontweight','b','FontSize',30)
% print('feedback_2lig_dust_export.png','-dpng')
if (model_run=='2L'),
    print('feedback_2L_dust_export.png','-dpng')
elseif (model_run=='2H'),
    print('feedback_2H_dust_export.png','-dpng')
end

figure(2)
h1 = plot(percentage*100, surf_fe_run1);
set(h1,'LineWidth',2)
hold on
h2 = plot(percentage*100, surf_fe_run2);
set(h2,'LineWidth',2)
set(gca,'FontSize',12)
xlabel('change in total dust (%)')
ylabel('global average surface dFe (nmol L^{-1})')
xl  = get(gca,'XLim');
yl  = get(gca,'YLim');
xpos = xl(1) + 0.85*(xl(2) - xl(1));
ypos = yl(1) + 0.1 *(yl(2) - yl(1));
ht = text(xpos,ypos,'b');
% set(ht,'Fontweight','b','FontSize',30)
print('feedback_2lig_dust_surf_fe.png','-dpng')
if (model_run=='2L'),
    print('feedback_2L_dust_surf_fe.png','-dpng')
elseif (model_run=='2H'),
    print('feedback_2H_dust_surf_fe.png','-dpng')
end

figure(3)
h1 = plot(percentage*100, av_fe_run1);
set(h1,'LineWidth',2)
hold on
h2 = plot(percentage*100, av_fe_run2);
set(h2,'LineWidth',2)
set(gca,'FontSize',12)
xlabel('change in total dust (%)')
ylabel('global average dFe (nmol L^{-1})')
xl  = get(gca,'XLim');
yl  = get(gca,'YLim');
xpos = xl(1) + 0.85*(xl(2) - xl(1));
ypos = yl(1) + 0.1 *(yl(2) - yl(1));
ht = text(xpos,ypos,'a');
set(ht,'Fontweight','b','FontSize',30)
% print('feedback_2lig_dust_ave_fe.png','-dpng')
if (model_run=='2L'),
    print('feedback_2L_dust_ave_fe.png','-dpng')
elseif (model_run=='2H'),
    print('feedback_2H_dust_ave_fe.png','-dpng')
end

figure(4)
h1 = plot(percentage*100, SO_fe_run1);
set(h1,'LineWidth',2)
hold on
h2 = plot(percentage*100, SO_fe_run2);
set(h2,'LineWidth',2)
set(gca,'FontSize',12)
xlabel('change in total dust (%)')
ylabel('Southern Ocean dFe (nmol L^{-1})')
xl  = get(gca,'XLim');
yl  = get(gca,'YLim');
xpos = xl(1) + 0.85*(xl(2) - xl(1));
ypos = yl(1) + 0.1 *(yl(2) - yl(1));
ht = text(xpos,ypos,'c');
set(ht,'Fontweight','b','FontSize',30)
%print('feedback_2lig_dust_SO_fe.png','-dpng')
if (model_run=='2L'),
    print('feedback_2L_dust_SO_fe.png','-dpng')
elseif (model_run=='2H'),
    print('feedback_2H_dust_SO_fe.png','-dpng')
end

