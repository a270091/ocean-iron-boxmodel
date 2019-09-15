function [f,po4_init,po4_final] = costf_mix_boxmodel_po4(pvec,rplot)
% f = costf_boxmodel_po4(pvec)

% initialize model parameters
global params
boxmodel_init_params()

% reminfrac(i,j) says which fraction of export from box j is remineralized
% in box i (all fractions are between zero and, and the sum of all
% fractions in a column is one)
params.reminfrac = zeros(12,5);
params.reminfrac(6,1) = 1 - pvec(1) - pvec(2);
params.reminfrac(8,1) = pvec(2);
params.reminfrac(11,1)= pvec(1);

params.reminfrac(6,2) = 1 - pvec(1) - pvec(2);
params.reminfrac(8,2) = pvec(2);
params.reminfrac(11,2)= pvec(1);

params.reminfrac(8,3) = 0.54;
params.reminfrac(9,3) = 0.06;
params.reminfrac(10,3)= 0.36;
params.reminfrac(11,3)= 0.02;
params.reminfrac(12,3)= 0.02;

params.reminfrac(7,4) = 0.85;
params.reminfrac(10,4)= 0.13;
params.reminfrac(12,4)= 0.02;

params.reminfrac(10,5)= 0.97;
params.reminfrac(12,5)= 0.03;

% initial PO4 distribution
po4_init = params.po4init;

% integrate
tspan = [0:50:3000];
po4 = ode23(@boxmodel_dgl_po4, tspan, po4_init);

% calculate RMS difference between final and initial PO4
po4_final = po4.y(:,end);
diff = po4_init - po4_final;
f = sqrt( sum( diff.^2 ) );

% prevent taking too strong mixing
if (max(abs(pvec))>40), 
    f = f + max(abs(pvec))-40;
end
if (min(pvec)<0), 
    f = f + 1000;
end

if (nargin>1),
    plot(po4.x,po4.y);
end

return
end
