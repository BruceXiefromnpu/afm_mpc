function [Y, U, U_nominal, dU, Xhat, Xerr] = sim_AFM(sim_struct, ref_traj) %#ok<STOUT>
% [Y, U, dU] = sim_MPC_fp(sim_struct, ref_f)
% 
% Inputs
% ------
%  sim_struct : a structure which must have the following feilds
%
%    (K_lqr | mpcProb);
%    PLANT;
%    du_max;
%    Nx_recyc;
%    sys_obs
% 
%  ref_traj : (timeseries) the setpoint to track
% Outputs
% -------
%  Y : (timeseries) Plant Output
%  U : (timeseries) Accumulated deltaU
%  dU : (timeseries) deltaU(k) Control (ie, output of MPC block)
%
% ------------------------------------------------------------------- %
% Pull out all the data stored in sim_struct to expose it to
% simulink. There must be a better way...

% tell the linter not to complain about unused variables, because they get used
% by simulink.
%#ok<*NASGU>
  VSS_LINEAR_CONTROL=Simulink.Variant('VSS_CONTROL_MODE==1');
  VSS_MPC_CONTROL=Simulink.Variant('VSS_CONTROL_MODE == 2');
  fprintf('-----------------------\n');  
  % Expose the sim struct to simulink.
  if isfield(sim_struct, 'mpcProb')
    mpcProb1 = sim_struct.mpcProb;
    VSS_CONTROL_MODE=2;
    % keyboard
    mpcProb1.reset_warm_start_data();
    fprintf('simulating as MPC\n')
  else
    VSS_CONTROL_MODE=1;
    fprintf('simulating as LINEAR\n');
  end
  PLANT = sim_struct.PLANT;
  trun = ref_traj.Time(end);
  Ts = PLANT.Ts;
  %trun = 5*Ts;
  x0 = SSTools.getNxNu(PLANT)*0;
  uss_0 = 0;
  
  
  
  sim_struct.x0_obs = sim_struct.sys_obs.b*0;
  [ ndist, Ns_obs] = size(sim_struct.sys_obs.c);
  Ident_obs = eye(Ns_obs);
    
  % the last row
  C_ydist = Ident_obs(end-ndist+1:end, :);
  % all rows but the last 1
  Ident_obs = Ident_obs(1:end-ndist, :);
  
  if isfield(sim_struct, 'thenoise')
    thenoise = sim_struct.thenoise; 
  else
    thenoise = timeseries(ref_traj.Time*0, ref_traj.Time);
  end
  
  if ~isfield(sim_struct, 'step_amp')
    sim_struct.step_amp = 0;
  end

  if isfield(sim_struct, 'r') && isfield(sim_struct, 'w')
    r = sim_struct.r;
    w = sim_struct.w;
  else 
    r = 0;
    w = 0;
  end
  
  if isfield(sim_struct, 'rp') && isfield(sim_struct, 'wp')
    rp = sim_struct.rp;
    wp = sim_struct.wp;
  else 
    rp = 0;
    wp = 0;
  end
     
  if ~isfield(sim_struct, 'gdrift_inv')
    sim_struct.gdrift_inv = tf(1,1, PLANT.Ts);
  end
  if ~isfield(sim_struct, 'gdrift')
    sim_struct.gdrift = tf(1,1, PLANT.Ts);
  end  
  
  % ------------------------  RUN THE SIM ---------------------------- %
  options = simset('SrcWorkspace','current');
  if ~sim_struct.isfxp
    fprintf('simulating as floating-point\n');
    sim('AFMss_fp_obshas_uk', [], options)
  else
    % Turn off warning about fxp states not getting logged.
    fprintf('simulating as FIXED-point\n');
    nw = sim_struct.nw;
    nf = sim_struct.nf;
    id = 'Simulink:Engine:BlkIgnoringUsedAsDStateFlag';
    warning('off', id);
    sim('FXP_AFMss_obshas_uk', [], options)
  end
  % provides Y, U, dU, Xhat, U_nominal
end
