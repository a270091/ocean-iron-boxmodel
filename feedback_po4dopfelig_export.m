%---------------------------------------------------------------------------------
% Feedback analysis for a box model with variable ligands
% 
% The faadback analysis requires three different runs of the box model:
% - the model with prognostic ligands is run into equilibrium for standard 
%   external Fe input and other model parameters
% - the equilibrium ligand concentration from the first run is extracted and
%   a second model run is done where the external forcing is varied, 
%   keeping the ligand distribution fixed at the distribution from the first run
% - the model is run a third time, with varied external forcing, but this time 
%   allowing the ligands to change with repect to the first run. 
%---------------------------------------------------------------------------------

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
conc_init = [po4_init;dop_init;fe_init;lig_init];

%---------------------------------------------------------------------------------
% run0: first integration of the box model, with variable ligands
%---------------------------------------------------------------------------------

tspan = (0:50:5000);
conc = ode23s(@boxmodel_dgl_po4dopfelig_export, tspan, conc_init);

% final state of the model, and some diagnostics
conc_final_run0 = conc.y(:,end);
c_export_run0 = export * params.redfield_c2p * 12 * 1.0e-18;

%---------------------------------------------------------------------------------
% changing parameters. Here it may be interesting to check a bit which ones one 
% would change. In this case we change total dust deposition
%---------------------------------------------------------------------------------

params.dust = params.dust*1.5; % increase total dust by XX percent

%---------------------------------------------------------------------------------
% run1: box model, with ligands kept fixed at distributions from run0,
% i.e. without the feedback
%---------------------------------------------------------------------------------

global ligfix
ligfix = conc_final_run0(37:48);

conc_init_nofback = conc_init(1:36);
tspan = (0:50:5000);
conc = ode23s(@boxmodel_dgl_po4dopfe_ligfix_export, tspan, conc_init_nofback);

% final state of the model, and some diagnostics
conc_final_run1 = conc.y(:,end);
c_export_run1 = export * params.redfield_c2p * 12 * 1.0e-18;

%---------------------------------------------------------------------------------
% run2: box model, with ligands reacting to the changed parameters, i.e. 
% with the ligand feedback included
%---------------------------------------------------------------------------------

tspan = (0:50:5000);
conc = ode23s(@boxmodel_dgl_po4dopfelig_export, tspan, conc_init);

% final state of the model, and some diagnostics
conc_final_run2 = conc.y(:,end);
c_export_run2 = export * params.redfield_c2p * 12 * 1.0e-18;

for k=1:12
  fprintf('Box %i: %s\n',k,params.long_names{k});
  fprintf('  DFe: %5.2f %5.2f %5.2f\n', conc_final_run0(k+24), conc_final_run1(k+24), conc_final_run2(k+24));
  fprintf('  PO4: %5.2f %5.2f %5.2f\n', conc_final_run0(k), conc_final_run1(k), conc_final_run2(k));
  fprintf('  DOP: %5.2f %5.2f %5.2f\n', conc_final_run0(k+12), conc_final_run1(k+12), conc_final_run2(k+12));
  fprintf('  Lig: %5.2f %5.2f %5.2f\n', conc_final_run0(k+36), ligfix(k), conc_final_run2(k+36));
  if k<=5,
    fprintf('  Export: %5.2f %5.2f %5.2f PgC \n', c_export_run0(k), c_export_run1(k), c_export_run2(k));
  end
  fprintf('\n');
end

