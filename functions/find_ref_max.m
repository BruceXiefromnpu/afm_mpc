% ref_max = find_ref_max(sim_struct, ref_s, varargin.)
%
% This function runs a (saturated) LQR or MPC simlation with the same
% control parameters over many different setpoints to determine the largest
% achievable setpoint (before instability results due to saturation).
% 
% It is designed to run the MPC_fp simulink model.
%
% Inputs
% ------
%   sim_struct: a struct containing the data simulink needs to run. It must
%               contain the following fields
%
%     sim_struct.PLANT :  the LTI plant model
%     sim_struct.Nx    :  normalized Xss
%     sim_struct.K_lqr :  LQR gain (equivilent gain if doing to MPC
%     simulation)
%     sim_struct.mpcProb1 : an instance of sparseMPCprob
%     sim_struct.trun : how long the simulation should for
%     sim_struct.du_max : the slew rate saturation limit
%     sim_struct.mpc_on : boolean flag of weather to run linear or MPC sim.
%
%     ref_s : a vector of references to check. Should be in
%     increasing order, as they will be tested from smallest to largest.
%
% Option Inputs
% ------------
%   find_ref_max(..., 'verbose', int) if int<=2, no plotting is down. If
%   int>2, will plot the y(t). If int > 3, will also plot for
%   u(t), and du(t) for each tested setpoint into the same set of
%   three figures. Default is verbose=0.
% 
%  find_ref_max(..., 'fig_base', int). Only used if verbose >2. Allows the
%  figures generated to be incremented. The scheme is 
%  figure(1 + fig_base) etc. Default is fig_base=10;
%
% Outputs
% -------
%   result.The maximim stable setpoint that was found. If the very first
%   iteration is unstable, ref_max=0.
%   result.ref_max : The maximum achievable reference.
% 
%   result.y_traj_s : A vector of time series trajectories from
%   ref_s(1) to ref_s(ref_max_ind).
% 
%   result.t_settle_s : A vector of all the settle times from
%   ref_s(1) to ref_s(ref_max_ind).
%
% See Also
% --------
%   sparseMPCprob, build_max_setpoints, StepData, StepParamsLin,
%   StepParamsMPC
%
% Depends on the simulink model 'MPC_fp'.

function result = find_ref_max(sim_struct, ref_s, varargin)
  % Auxiliary functions for checking parameter values:
  isnonnegint = @(x) (isnumeric(x) && mod(x, 1) == 0 && x >= 0);
  
  % Parse options
  defaultVerbose = 0;
  defaultFigbase = 10;
  p = inputParser;
  addParameter(p, 'verbose', defaultVerbose, isnonnegint)
  addParameter(p, 'fig_base', defaultFigbase, isnonnegint)
  addParameter(p, 'max_sp_judge', @max_sp_judge_default);
  parse(p, varargin{:})
  verbose = p.Results.verbose;
  fig_base = p.Results.fig_base;
  max_sp_judge = p.Results.max_sp_judge;
  
  figs = [];
  if verbose >= 2
    F1 = figure(1 + fig_base);
    figs = [figs, F1];
  end
  if verbose >= 3
    F2 = figure(2 + fig_base);
    F3 = figure(3 + fig_base);
    figs = [figs, F2, F3];
  end
  
  % We have to set the simulink workspace to this function, since
  % default is the main scripts workspace.
  % options = simset('SrcWorkspace','current');
  y_traj_s = repmat(timeseries(0,0), 1, length(ref_s));
  du_traj_s = repmat(timeseries(0,0), 1, length(ref_s));
  t_settle_s = zeros(1, length(ref_s));
  ref_max_recommended_found = 0;
  ref_max_recommended_idx = 0;
  for iter = 1:length(ref_s)
    ref_f = ref_s(iter);
    % sim('MPC_fp', [], options)
    % y1 = y_mpcDist;
    [Y, U, dU] = sim_MPC_fp(sim_struct, ref_f);
    [t_settle, k_s] = settle_time(Y.time, Y.Data, ref_f, 0.01*ref_f,...
      [], [], 30);
    
    if verbose
      figs = plot_local(Y, U, dU, figs, fig_base, verbose);
    end
    
    if ~ref_max_recommended_found
      ref_max_recommended_found = max_sp_judge(ref_f, ...
        iter, t_settle, Y); % is =1 if we violate threshold, zero otherwise
      % If we don't violate threshold, add one to the idx counter. Otherwise,
      % add nothing. 
      ref_max_recommended_idx = ref_max_recommended_idx + ~ref_max_recommended_found;
      %             ref_max_recommended_idx = max(1, iter-1);
    end
    if isnan(t_settle)
      ts_is_nan = true;
      break
    else
      ts_is_nan = false;
      t_settle_s(iter) = t_settle;
      y_traj_s(iter) = Y;
      du_traj_s(iter) = dU;
    end
  end
  
  % if ~ref_max_recommended_found
  %     ref_max_recommended_idx = iter;
  % end
  
  
  if ts_is_nan && iter == 1 % The smallest reference was unstable.
    result.ref_max_idx = NaN;
    result.ref_max_recommended_idx = NaN;
    result.t_settle_s = [];
    result.y_traj_s = Y;
    result.du_traj_s = dU;
    return
  elseif ts_is_nan % we reached an unstable ref, but its not the
    % first. So largest stable is the last one.
    ref_max_idx = iter-1;
  else             % We never found an unstable ref.
    ref_max_idx = iter;
  end
  result.ref_max_idx = ref_max_idx;
  result.ref_max_recommended_idx = ref_max_recommended_idx;
  result.y_traj_s = y_traj_s;
  result.du_traj_s = du_traj_s;
  result.t_settle_s = t_settle_s;
  
end

function ref_max_recommended_found = max_sp_judge_default(varargin)
  
  ref_max_recommended_found = false;
end


function figs = plot_local(y1, u1, du1, figs, fig_base, verbose)
  % Figure (1)
  if verbose >= 2
    if isvalid(figs(1))
      change_current_figure(figs(1))
    else
      figs(1) = figure(1 + fig_base);
    end
    hold on
    plot(y1.Time, y1.Data)
    title('y(t)')
  end
  
  if verbose >= 3
    % Figure 2
    if isvalid(figs(2))
      change_current_figure(figs(2))
    else
      figs(2) = figure(2 + fig_base);
    end
    hold on
    plot(u1.Time, u1.Data)
    title('u(t)')
    
    % Figure 3
    if isvalid(figs(3))
      change_current_figure(figs(3))
    else
      figs(3) = figure(3 + fig_base);
    end
    hold on
    plot(du1.Time, du1.Data)
    title('du(t)')
  end
  drawnow()
  
  
end








