clear, clc

Ts = 40e-6;
dry_run = true;
saveon = true;

step_amp = 0.15;
step_start_idx = 3000;
step_end_idx = 10500;
u = zeros(step_end_idx, 1);
u(step_start_idx) = step_amp;
u_vec = cumsum(u);

t_vec = (0:length(u_vec)-1)'*Ts;
figure(1); clf
plot(t_vec, u_vec);
grid on



reset_piezo('t1', 15, 't_final', 25, 'umax', 3, 'k1', 0.2,...
            'verbose', true, 'dry_run', dry_run)
if ~dry_run
  clear vi;
  vipath_reset = 'C:\Users\arnold\Documents\MATLAB\afm_mpc_journal\labview\reset_piezo.vi';
      [e, vi] = setupVI(vipath_reset, 'Abort', 0,...
    'umax', 9, 'TsTicks', 1600, 'u_in', u_vec);
  vi.Run;
  stage_dat = vi.GetControlValue('stage_data_out');

  u_exp = stage_dat(:,1);
  y_exp = stage_dat(:,2); % - dat(1,2);
  t_exp = (0:length(u_exp)-1)'*Ts;

  figure(201); clf
  hold on, grid on
  plot(t_exp, u_exp)
  plot(t_exp, y_exp - y_exp(1))
  
  if saveon
    driftData.t_exp = t_exp;
    driftData.u_exp = u_exp;
    driftData.y_exp = y_exp;
    save('driftID_data_4-30-2018_01.mat', 'driftData')
  end
end


ms = 1e3;








