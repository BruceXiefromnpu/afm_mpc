


clear, clc
addpath('many_steps_data')
load('many_steps.mat')

whos
L = 800;
TOL = 0.02;

ref_s = ref_traj_params.ref_s;
step_idx = ref_traj_params.impulse_idx;

load('many_steps_hyst_withsat.mat')
TS_hystsat = get_many_steps_ts(y_exp, ref_s, step_idx, TOL, 1);

figure(1000)
plot(ref_traj_params.ref_traj.Time, ref_traj_params.ref_traj.Data)
hold on
h1 = plot(y_exp.Time, y_exp.Data, '-r');
h1.DisplayName = 'linfp, inv hyst w/ sat ';
%

load('many_steps_hyst_nosat_R4.mat')
TS_hyst = get_many_steps_ts(y_exp, ref_s, step_idx, TOL, 1);
h2 = plot(y_exp.Time, y_exp.Data, '-b');
h2.DisplayName = 'Invert hyst (linfp, no sat operator)';
%
load('many_steps_noinvert.mat')
TS_noinv = get_many_steps_ts(y_exp, ref_s, step_idx, TOL, 1);
h3 = plot(y_exp.Time, y_exp.Data, '-g');
h3.DisplayName = 'linfp, no inversion';
%
load('many_steps_pi.mat')
TS_pi = get_many_steps_ts(y_exp, ref_s, step_idx, TOL, 1);
h4 = plot(y_exp.Time, y_exp.Data, '-m');
h4.DisplayName = 'PI';


load('many_steps_mpc_invHyst_invDrift.mat')
TS_mpc = get_many_steps_ts(y_exp, ref_s, step_idx, TOL, 1);
h5 = plot(y_exp.Time, y_exp.Data, '--k');
h5.DisplayName = 'MPC (inv H \& drft)';
leg = legend([h1, h2, h3,h4, h5]);
set(leg, 'Location', 'NorthEast')

%%

S = sprintf('\\begin{tabular}{cccccccc}\n')
S = sprintf(['%s&ref & delta & Ts (no inv) & Ts (inv H+S) ',...
             '& Ts (inv H) & Ts (PI) & TS (MPC)\\\\\n'], S);
S = sprintf('%s\\toprule\n', S);

for k = 2:length(ref_s)
  
  delta_ref = ref_s(k) - ref_s(k-1);
  S = sprintf('%s &%.2f & %.2f & %.3f & %.3f & %.3f & %.3f & %.3f\\\\\n', S, ref_s(k), delta_ref,...
    1000*TS_noinv(k-1), 1000*TS_hystsat(k-1), 1000*TS_hyst(k-1),...
    1000*TS_pi(k-1), 1000*TS_mpc(k-1));
end

S = sprintf('%s\\midrule\n total & -- & --& %.2f & %.2f &%.2f & %.2f & %.2f\\\\\n', S,...
               1000*sum(TS_noinv), 1000*sum(TS_hystsat), 1000*sum(TS_hyst),...
               1000*sum(TS_pi), 1000*sum(TS_mpc))

S = sprintf('%s\\end{tabular}\n', S);
%%
fid = fopen('Z:\manystepsdata.tex', 'w+')

fprintf(fid, '%s', S)

fclose(fid)
