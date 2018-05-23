clear, clc

% 1. v_low delta max. From G_deltauz2Ipow
load(fullfile(PATHS.sysid, 'FRF_data_current_stage2.mat'));

G_deluz2Ipow = modelFit.models.G_deluz2powI;
I_max = StageParams.Imax;

% compute impulse response for 1-norm:
[y, t] = impulse(G_deluz2Ipow);
y = y*G_deluz2Ipow.Ts; % matlab makes the impulse 1/Ts tall.
nm1 = sum(abs(y))

du_max = I_max/nm1



%%

models = CanonPlants.plants_with_drift_inv(true)

gd_inv = models.gdrift_inv;

[y, t] = impulse(gd_inv);
y = y*gd_inv.Ts;

nm1_drift = sum(abs(y))

% compare to the retardation I was doing before
norm(gd_inv)

du_max_effective = du_max/nm1_drift

