% Build the u-reset.
%
% Options
% ----------
%
% t1, t_final, umax, k1, verbose, dry_run
function [stage_data, u_reset] = reset_piezo(varargin)
P = path();
  addpath('C:\Users\arnold\Documents\MATLAB\afm_mpc_journal\modelFitting\hysteresis')

  p = inputParser;
  p.addParameter('t1', 10)
  p.addParameter('t_final', 20);
  p.addParameter('umax', 5);
  p.addParameter('k1', 0.45);
  p.addParameter('verbose', true);
  p.addParameter('dry_run', false);
  
  p.parse(varargin{:});
  t1 = p.Results.t1;
  t_final = p.Results.t_final;
  umax = p.Results.umax;
  k1 = p.Results.k1;
  verbose = p.Results.verbose;
  dry_run = p.Results.dry_run;
  Ts = 40e-6;
  % t1 = 10;
  % t_final = 20;
  % umax = 5;
  % k1 = 0.45;
  
  u_reset = PIHyst.gen_reset_u(t1, t_final, Ts, k1, umax);
  if verbose
    figure(100);
    t = (0:length(u_reset)-1)'*Ts;
    plot(t, u_reset, '-k')
    grid on, hold on
  end
  slewfname_in = 'hyst_reset_datain.csv';
  slewfpath_in = fullfile(pwd, slewfname_in);
  slewfname_out = 'hyst_reset_data_out.csv';
  slewfpath_out = fullfile(pwd, slewfname_out);
  
  
  if ~dry_run
    clear vi;
    vipath_reset = 'C:\Users\arnold\Documents\MATLAB\afm_mpc_journal\labview\reset_piezo.vi';
        [e, vi] = setupVI(vipath_reset, 'Abort', 0,...
      'umax', umax+0.25, 'TsTicks', 1600, 'u_in', u_reset);
    vi.Run;
    stage_dat = vi.GetControlValue('stage_data_out');

    if verbose
      u_exp = stage_dat(:,1);
      yx_exp = stage_dat(:,2); % - dat(1,2);
      t_exp = (0:length(u_exp)-1)'*Ts;

      figure(100); 
      plot(t_exp, u_exp, 'r')
      hold on
      plot(t_exp, yx_exp, 'b')
      grid on
    end
  end
  path(P)
end
