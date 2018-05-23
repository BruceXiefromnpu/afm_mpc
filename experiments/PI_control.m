clear
load(fullfile(PATHS.sysid, 'FRF_data_current_stage2.mat'))
save_root = fullfile(PATHS.exp, 'step-exps', 'many_steps_data_rand');
TOL = 0.01;

G = ss(modelFit.models.G_uz2stage);
G = absorbDelay(G);
Ts = G.Ts;
z= tf('z', Ts);
% Ki = 0.02;
% Kp = 0.05;
% D = zpk(Kp + Ki/(z-1))
C = 0.035972*(z-0.5298)/(z-1)
D=C;
[num,den] = tfdata(D, 'v');
Kp = num(1);
Ki = num(2) + Kp;

H = feedback(G*D, 1);
H2 = feedback(G*C, 1);
% figure(22); clf
% t = (0:2000)'*Ts;
% step(H,H2, t)
% grid on

N1    = 800;
ref_f_1 = 1; % 1.5 to hit slew rate, 1.4 doesn't
trajstyle = 3;
if trajstyle == 1
    N2 = 1200;
    trun = Ts*(N1 + N2);
    ref_f_2 = -4; % 1.5 to hit slew rate, 1.4 doesn't
    ref_0 = 0;
    t1 = [0:1:N1]'*Ts;
    t2 = [N1+1:1:N1+N2]'*Ts;

    t = [t1; t2];
    yref = [0*t1 + ref_f_1;
            0*t2 + ref_f_2];
    ref_traj = timeseries(yref, t);
elseif trajstyle == 2
    t1 = [0:1:N1]'*Ts;
    trun = Ts*(N1);
    t = t1;
    yref = [0*t1 + ref_f_1];
    ref_0 = 0;
    ref_traj = timeseries(yref, t);
elseif trajstyle == 3
  
  load(fullfile(save_root, 'many_steps_rand_longts.mat'))
  ref_traj = ref_traj_params.ref_traj;
end


PLANT = G;
x0 = G.b*0;
% D=C;
trun = ref_traj.Time(end)
sim('AFMss_TFcontrol')

piOpts = stepExpOpts('pstyle', '-r', 'TOL', TOL, 'y_ref', ref_f_1,...
                      'controller', D, 'name',  'FP Simulation');
                    
                    
sim_pi = stepExpDu(Y, U, dU, piOpts);

F1 = figure(60); clf
plot(sim_pi, F1, 'umode', 'both');


ref_s = ref_traj_params.ref_s;
step_idx = ref_traj_params.impulse_idx;
TS_pi = get_many_steps_ts(Y, ref_s, step_idx, TOL, 1);

%%
traj_path = 'C:\Users\arnold\Documents\MATLAB\afm_mpc_journal\data\traj_data.csv';
dat_out_path = 'Z:\mpc-journal\step-exps\expOut_PI.csv';
fid = fopen(traj_path, 'w+')
fprintf(fid, '%.6f, ', ref_traj.Data);
fprintf(fid, '\n');
fclose(fid);
Iters = length(ref_traj.Data);
umax = 10;
ymax = 8;

SettleTicks = 1000;
% ymax = max(abs(ref_traj.Data(:)))*1.2;

Iters = min(Iters, length(ref_traj.data));

vipath ='C:\Users\arnold\Documents\MATLAB\afm_mpc_journal\labview\play_AFM_PI.vi';

[e, vi] = setupVI(vipath, 'SettleTicks', SettleTicks, 'Iters', Iters,...
            'umax', umax, 'ymax', ymax, 'outputDataPath', dat_out_path,...
            'traj_path', traj_path, 'Ki', Ki, 'Kp', Kp);
          
vi.Run
%
AFMdat = csvread(dat_out_path);
t_exp = (0:size(AFMdat,1)-1)'*Ts;
y_exp = timeseries(AFMdat(:,1), t_exp);
u_exp = timeseries(AFMdat(:,2), t_exp);
I_exp = timeseries(AFMdat(:,4), t_exp);

du_ = diff(u_exp.Data);
du_exp = timeseries([du_; du_(end)], t_exp);

pi_expOpts = stepExpOpts('pstyle', '-g', 'TOL', TOL, 'y_ref', ref_f_1,...
                      'controller', D, 'name',  'FP Simulation');
                    
                    
sim_exp_pi = stepExpDu(y_exp, u_exp, du_exp, pi_expOpts);

plot(sim_exp_pi, F1, 'umode', 'both'); 


ref_s = ref_traj_params.ref_s;
step_idx = ref_traj_params.impulse_idx;
TS_pi = get_many_steps_ts(y_exp, ref_s, step_idx, TOL, 1, 'abs');
%%
save(fullfile(save_root, 'many_steps_pi.mat'), 'y_exp', 'u_exp', 'I_exp', 'ref_s', 'ref_traj')

