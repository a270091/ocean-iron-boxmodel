function [f,po4_init,po4_final] = costf_remin_boxmodel_po4dopfe_export(pvec,rplot)

% initialize model parameters
global params
boxmodel_init_params()

if (nargin>0) 
  remfrac = params.reminfrac;
  remfrac(8,3) = pvec(1);
  remfrac(9,3) = 0.2;
  remfrac(10,3)= 1 - pvec(1) - 0.24;
  remfrac(11,3)= 0.02;
  remfrac(12,3)= 0.02;
  params.reminfrac = remfrac;
end

% initial PO4 distribution
po4_init = params.po4init;
dop_init = zeros(size(po4_init));
fe_init = 0.6 + zeros(size(po4_init));
conc_init = [po4_init;dop_init;fe_init];

% integrate
tspan = (0:50:3000);
conc = ode23s(@boxmodel_dgl_po4dopfe_export, tspan, conc_init);

% calculate RMS difference between final and initial PO4
po4_final = conc.y(1:12,end);
diff = po4_init - po4_final;
% f = sqrt( sum( (diff.^2).*params.volume ) ./ sum(params.volume) );
f = sqrt( sum( (diff.^2) ) );
%f = sum( abs(diff ) );

if (nargin>1)
    if (rplot>0)
        plot(conc.x,conc.y);
    end
end

return
end
