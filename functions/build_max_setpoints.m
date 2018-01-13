
function step_data = build_max_setpoints(step_data, varargin)
% step_data = build_max_setpoints(step_data, varargin) 
%    
% Builds a set of maximum setpoints for a family of LQR based
% linear or MPC based controllers. The controllers are
% parameterized by gamma. 
%
% Inputs
% -----
%   step_data : a class instance of StepData 
%   
% Outputs
% -------
%   step_data : the function will populate the StepData.results
%   property upon exit. In addition, if the StepData.file field has
%   been set, will also save the resulting instance into that file
%   location.
% 
%   Upon (succseful exit), the results property will contain the
%   following fields for the linear case:
%   results.max_setpoints.
%          results.data{idx}.t_settle_s
%          results.data{idx}.y_traje_s
%          results.data{idx}.ref_max
%   where idx is the index for gamma.
% 
%   For the mpc case, results is itself a cell array, with size=length(N_mpc).
% 
% Optional Inputs
% --------------
%   build_max_setpoints(..., 'force', (true|false)) if true,  will
%   force the simulations to be re-run regardless of if the
%   underlying parameters are the same as the saved ones.
% 
%   build_max_setpoints(..., 'fid', fid) Default fid=1. Pass in a
%   file id from fopen to write all logging info to a file. Will
%   only be used if step_data.verbose >0.
%
%   verbose = 0 --> only initial and final status printed to console.
%   verbose =  1 print progress bar to console
%   verbose =  2 plot gamma vs max in Fig 1, and settle-times vs
%   ref in Fig (all gammas).
%   verbose = 3 plot y-trajs into a figure. A new figure is
%   generated for each gamma.
%   verbose = 4 also plot u/du-trajs into 2 figures. A new figure is
%   generated for each gamma. This generates a lot of figures, so
%   currently, they are all deleted after every 10th gamma.
% 
% Algorithm/Methods
% ----------------
%   Linear + delU Saturation Case 
%   ----------------------------- 
%    Iterate over a bunch of gammas and try to find the maximum setpoint for
%    each one. This code is pretty niave. We start at a very low setpoint and
%    slowly increase the setpoint until the settling time is reported as NaN.
%    A bisection search would probably be faster. This is complicated though
%    by the fact that instability seems to occur in the "middle": if the
%    setpoint is large enough, we dont have the stability problem, which is
%    weird.
%
%  MPC Case 
%  --------
%    This works similarly, but now we also allow the simulations to
%    be done over a family of control horizons. Thus, there are two
%    loops: the outer loop loops of N_mpc_s, the inner loop loops
%    over gammas, like in the linear case.
    
    
    % pull data out of data variable:
    defaultForce = 0;
    p = inputParser;
    p.addParameter('force', defaultForce);
    p.addParameter('fid', 1);
    p.addParameter('verbose', 1);
    parse(p, varargin{:});
    force = p.Results.force;
    fid = p.Results.fid;
    verbose = p.Results.verbose;
    if stepdata_struct_unchanged(step_data) && ~force
        load(step_data.file)
        fprintf(fid, 'LOG (build_max_setpoints, mpc_on=%d)\n', step_data.params.sim_struct.mpc_on);
        fprintf(fid, ['Data appears to be the same. Loading data ',...
                 'without re-calculation.\n']);
        return
    end
    fprintf(fid, 'LOG (build_max_setpoints, mpc_on=%d)\n', step_data.params.sim_struct.mpc_on);
    fprintf(fid, 'data has changed: re-building max setpoints.\n');
    
    params = step_data.params;
    % Expose parameters:
    
     
    gam_s = params.gam_s;
    ref_s = params.ref_s;
    du_max = params.du_max;
    
    mpc_on = params.sim_struct.mpc_on;

    if mpc_on
        max_sp_fieldname = 'max_setpoints_mpc';
        N_mpc_s = params.N_mpc_s;
        mpc_mode = params.mpc_mode;
    else
        max_sp_fieldname = 'max_setpoints_lin';
    end
    
    if verbose > 1
        F1 = figure(1); clf;
        xlabel('$\gamma$', 'interpreter', 'latex', 'FontSize', 16)
        ylabel('Max Ref', 'interpreter', 'latex', 'FontSize', 16)
        
        F2 = figure(2); clf; hold on;
        ylabel('settle time [ms]')
        xlabel('ref_f')
        Figs = [F1, F2];
    else
        Figs = [];
    end
    
    % ------------------------------------------------------------------ %
    % Finally, run the simulations. Pass all this shit into the local
    % runner function.
    PB = [];
    
    if ~mpc_on
        if verbose > 0
            PB = ProgressBar(length(params.gam_s), 'start_str', ...
                             'Lin', 'fid', fid);
        end
        % Make it a cell array of size one for compatibility with
        % the MPC case.
        result_s{1} = build_max_sp_local(params, Figs, ...
                                           PB, verbose);
    else
        %      result_s.max_setpoints = max_setpoints;
        %      result_s.data = result_data; 
        %        -Where max_setpoints is a vector (list) of the maximum
        %         setpoint achievable for each supplied gamma.
        %        - data is a vector of structs containing:
        %            data(idx).t_settle_s
        %            data(idx).y_traj_s
        %            data(idx).ref_max
        
        result_s = cell(1, length(N_mpc_s));
        sim_struct.mpc_on = true;
        for mpc_iter = 1:length(params.N_mpc_s)
            N_mpc_iter = params.N_mpc_s(mpc_iter);
            if verbose > 0
                start_str = sprintf('MPC, N=%.0f', N_mpc_iter);
                PB = ProgressBar(length(params.gam_s), 'start_str', ...
                                 start_str, 'fid', fid);
            end
 
            result_s_iter = build_max_sp_local(params, Figs, PB, verbose, ...
                                                           N_mpc_iter);
            result_s{mpc_iter} = result_s_iter;
        end
    
    end

    % ---------------- Save Data and Pretty-fy Figures---------------------%
    
    % results.max_setpoints = max_setpoints;
    step_data.results = result_s;
        

    if step_data.savedata
        save(step_data.file, 'step_data')
        if verbose >1
            saveas(F1, step_data.fig_files(1))
            saveas(F1, step_data.fig_files(2))
        end
    end

end % END MAIN FUNCTION

% -----------------------------------------------------------------------%
% -----------------------------------------------------------------------%

function result_s = build_max_sp_local(params, Figs, PB, verbose, ...
                                                 varargin)
% Returns
% ------
%      result_s.max_setpoints = max_setpoints;
%      result_s.data = result_data; 
%        -Where max_setpoints is a vector (list) of the maximum
%         setpoint achievable for each supplied gamma.
%        - data is a vector of structs containing:
%            data{idx}.t_settle_s
%            data{idx}.y_traj_s
%            data{idx}.ref_max
%         In other words, the structure contained in each element
%         of data contains the max reference achievable for
%         gamma(idx), and the y_trajectories and t_settle_s contain
%         the eponymous data for each reference tested up until
%         (and including) ref_max.
    
    sim_struct = params.sim_struct;
    if length(varargin) == 1
        N_mpc_iter = varargin{1};
    else
        N_mpc_iter = [];
    end
    

    gam_s = params.gam_s;
    ref_s = params.ref_s;
% Initialize the "warm starting"
    ref_max = ref_s(1);
    max_setpoints_iter = 0*gam_s;
    if verbose
        PB.upd(0);
    end
    data = cell(1, length(gam_s));
    for gam_iter = 1:length(gam_s)
        gamma = gam_s(gam_iter);
        % Update either K or the mpcProb.
        sim_struct = update_sim_struct(sim_struct, params, gamma, N_mpc_iter);
        
        % Warm starting
        rmax_ind = find(ref_s == ref_max, 1);
        rmax_ind = max([rmax_ind-5, 1]);
        
        fig_base = 10*gam_iter;
        data_iter = find_ref_max(sim_struct, ref_s(rmax_ind:end),...
                                             'verbose', verbose,...
                                             'fig_base', fig_base);
        max_setpoints(gam_iter) = data_iter.ref_max;
        data{gam_iter} = data_iter;
        
        if verbose > 1
            plot_local(data_iter.t_settle_s, gam_s, ref_s, max_setpoints, ...
                       Figs, gam_iter)
        end
        if verbose > 0
            PB.upd(gam_iter);
        end
        
        %     keyboard
    end % end MAIN LOOP 
    result_s.max_setpoints = max_setpoints;
    result_s.data = data;
end

function plot_local(t_settle_s, gam_s, ref_s, max_setpoints, Figs, gam_iter)
% expose for easy access while plotting:
    
    if mod(gam_iter, 10)==0
        fprintf('\ndeleting current figures\n\n')
        for i=1:gam_iter
            close(figure(10*i + 1))
            close(figure(10*i + 2))
            close(figure(10*i + 3))
        end
    end
    change_current_figure(Figs(2));
    % colrs = get(gca, 'colororder');
    k_max = find(t_settle_s ~= 0, 1, 'last');
    
    plot(ref_s(1:k_max), t_settle_s(1:k_max)*1e3)
    % ylim([0, 10])
    hold on
    drawnow()

    change_current_figure(Figs(1))
    hlin = plot(gam_s(1:gam_iter), max_setpoints(1:gam_iter),...
                '-o',  'LineWidth', 2);

    drawnow()
    hold on
end


function sim_struct = update_sim_struct(sim_struct, params, gamma, N_mpc_iter)
   if isempty(N_mpc_iter) % we are in linear mode then.
    sim_struct.K_lqr = dlqr(params.sys.a, params.sys.b,...
                                params.Q, gamma);
   else % if N_MPC_iter not empty --> in MPC mode.
       du_max = sim_struct.du_max;
       N_mpc  = N_mpc_iter;
       Qp = dare(params.sys.a, params.sys.b, params.Q, gamma); 
       sim_struct.mpcProb1 = condensedMPCprob(params.sys, N_mpc,...
                                              params.Q, Qp, gamma);
       sim_struct.mpcProb1.add_U_constraint('box', [-du_max, du_max]);
       sim_struct.K_lqr = params.sys.C*0;
   end
end



function status = stepdata_struct_unchanged(data_struct)
% This function will return 1 if
% (a) data_struct.file exists
% (b) the params fields in the loaded data struct (at data_struct.file) 
% and the provided data are the same. 
    
    if ~exist(data_struct.file, 'file')
        status = 0;
        % keyboard
        return 
    else
        load(data_struct.file)
        % Provides: step_data
    end
    params_other = step_data.params;
    if ~isequal(data_struct.params, params_other)
        status = 0;
        return
    end
    
    status = 1;
end


