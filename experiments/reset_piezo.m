% Build the u-reset.
function reset_piezo(varargin)
P = path();
  addpath('C:\Users\arnold\Documents\MATLAB\afm_mpc_journal\modelFitting\hysteresis')

  p = inputParser;
  p.addParameter('t1', 10)
  p.addParameter('t_final', 20);
  p.addParameter('umax', 5);
  p.addParameter('k1', 0.45);
  p.addParameter('verbose', true);
  
  p.parse(varargin{:});
  t1 = p.Results.t1;
  t_final = p.Results.t_final;
  umax = p.Results.umax;
  k1 = p.Results.k1;
  verbose = p.Results.verbose;
  Ts = 40e-6;
  % t1 = 10;
  % t_final = 20;
  % umax = 5;
  % k1 = 0.45;
  
  u_reset = PIHyst.gen_reset_u(t1, t_final, Ts, k1, umax);
  if verbose
    figure(100)
    plot(u_reset)
    grid on
  end
  slewfname_in = 'hyst_reset_datain.csv';
  slewfpath_in = fullfile(pwd, slewfname_in);
  slewfname_out = 'hyst_reset_data_out.csv';
  slewfpath_out = fullfile(pwd, slewfname_out);
  
  csvwrite(slewfpath_in, u_reset);
  clear vi;
  vipath_reset = 'C:\Users\arnold\Documents\MATLAB\afm_mpc_journal\labview\play_nPoint_id_slew_OL_FIFO.vi'
  
  [e, vi] = setupVI(vipath_reset, 'Abort', 0,...
    'umax', umax+0.25, 'data_out_path', slewfpath_out,...
    'traj_in_path', slewfpath_in, 'TsTicks', 1600);
  vi.Run
  if verbose
    dat = csvread(slewfpath_out);
    
    u_exp = dat(:,1);
    yx_exp = dat(:,2); % - dat(1,2);
    t_exp = (0:length(u_exp)-1)'*Ts;
    
    figure(101);
    plot(t_exp, u_exp*.68)
    hold on
    plot(t_exp, yx_exp)
    grid on
  end

  path(P)
end
