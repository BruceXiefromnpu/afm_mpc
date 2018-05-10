


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
h1.DisplayName = 'Invert hyst with sat operator';
%%

load('many_steps_hyst_nosat.mat')
TS_hyst = get_many_steps_ts(y_exp, ref_s, step_idx, TOL, 1);
h2 = plot(y_exp.Time, y_exp.Data, '-b');
h2.DisplayName = 'Invert hyst (no sat operator)';
%%
load('many_steps_noinvert.mat')
TS_noinv = get_many_steps_ts(y_exp, ref_s, step_idx, TOL, 1);
h3 = plot(y_exp.Time, y_exp.Data, '-g');
h3.DisplayName = 'no inversion';
%%

S = sprintf('\\begin{tabular}{cccccc}\n')
S = sprintf('%s&ref & delta & Ts (no inv) & Ts (inv H+S) & Ts (inv H)\\\\\n', S);
S = sprintf('%s\\toprule\n', S);

for k = 2:length(ref_s)
  
  delta_ref = ref_s(k) - ref_s(k-1);
  S = sprintf('%s &%.2f & %.2f & %.3f & %.3f & %.3f\\\\\n', S, ref_s(k), delta_ref,...
    1000*TS_noinv(k-1), 1000*TS_hystsat(k-1), 1000*TS_hyst(k-1));
end

S = sprintf('%s\\midrule\n total & -- & --& %.2f & %.2f &%.2f\\\\\n', S,...
               1000*sum(TS_noinv), 1000*sum(TS_hystsat), 1000*sum(TS_hyst))

S = sprintf('%s\\end{tabular}\n', S);
%
fid = fopen('Z:\manystepsdata.tex', 'w+')

fprintf(fid, '%s', S)

fclose(fid)
