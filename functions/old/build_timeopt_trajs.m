% clqr_data = build_clqr_trajs(data_struct)
%
% Builds CLQR trajectories parameteriezed by gamma. Will run a sequence
% of simulations over a list of references, for each gamma. Will record
% both the resulting trajectories as well as the associated settling
% times. 
%
% Inputs
% -----
%    data_struct : is a structure with the followin format:
%    data_struct.sys = sys_recyc;
%    data_struct.Q   = Q;
%    data_struct.ref_s = ref_s;
%    data_struct.sim_struct = sim_struct;
%    data_struct.verbose = 1;
%    data_struct.savedata = 1;
%    data_struct_clqr.gam_s_clqr = [200]; % The list of gammas to loop over.
%    data_struct_clqr.file = 'data/clqr_data.mat'; % The file to save
%                                                  % things to
%    data_struct_clqr.du_max = du_max : slew rate limit
%    data_struct_clqr.N_traj = 400    : control horizon.
%
% The big function of this function is to check if data_struct.file
% exists. If it does, then the function compares the data saved in that
% .mat file to the data provided in data_struct. If the data is the same,
% the function simply loads the .mat file and returns the structure
% contained therin. If the data has changed, the function will perform
% the simulations again, saving the new data into the .mat file specified
% in data_struct.file.

function [step_data, status] = build_timeopt_trajs(step_data, varargin)
% Description.
    defaultForce = 0;
    p = inputParser;
    p.addParameter('force', defaultForce);
    p.addParameter('max_iter', 20);
    p.addParameter('do_eject', true);
    p.addParameter('fid', 1);
    p.addParameter('verbose', 1);
    
    parse(p, varargin{:});
    force = p.Results.force;
    max_iter = p.Results.max_iter;
    do_eject = p.Results.do_eject;
    fid = p.Results.fid;
    verbose = p.Results.verbose;
    
    status = 0;
    if stepdata_struct_unchanged(step_data) && ~force
        load(step_data.file)
        fprintf(fid, 'LOG: (build_timeopt_trajs)\n');
        fprintf(fid, ['Data appears to be the same. Loading data ',...
                 'without re-calculation.\n\n']);
        return
    end
    fprintf(fid, 'LOG (build_timeopt_trajs:\n)');
    fprintf(fid, 'data has changed: re-building max setpoints.\n');
    
    params = step_data.params;
    
    % Expose parameters:
    ref_s = params.ref_s;
    du_max = params.du_max;
       
    sys_nodelay = params.sys_nodelay;
    Nd = params.Nd;
    
        
    if verbose
        Fig = figure(200);
        hands = [];
        change_current_figure(Fig);
        colrs = get(gca, 'colororder');
    end

    if do_eject
        % Pull out open-loop pole-zero information.
        Ts = sys_nodelay.Ts;
        [wp_real_x, wz_real_x] = w_zp_real(params.sys_nodelay);
        rho_1 = wz_real_x(1)/wp_real_x(1);
        g_eject = zpk(exp(-wp_real_x(1)*Ts), exp(-wz_real_x(1)*Ts), 1, Ts);
        g_eject = g_eject/dcgain(g_eject);
        sys_eject = minreal(sys_nodelay*g_eject);
        sys_sim = SSTools.deltaUkSys(sys_eject);
    else
        sys_sim = SSTools.deltaUkSys(sys_nodelay);
    end
    
    Nx_sim = SSTools.getNxNu(sys_sim);
    x0_sim = Nx_sim*0;
    
    % Pre-allocate
    
    opt_trajs_save = cell(1, length(ref_s));
    time_opt_settletime_s = ref_s*0;
    
    warning('off', 'MATLAB:nargchk:deprecated')
    toBisect = TimeOptBisect(sys_sim, du_max);
    toBisect.max_iter = max_iter;
    
    for iter = 1:length(ref_s)
        fprintf(fid, ['Time Optimal bisection for ref=%.3f, iter = %.0f ' ...
                 'of %.0f\n'], ref_s(iter), iter, length(ref_s));
        xf = Nx_sim*ref_s(iter);
        [X, U, status]=toBisect.time_opt_bisect(x0_sim, xf);
        if status
            fprintf(fid, ['Time-optimal bisection failed at ref = ' ...
            '%0.3f. Exiting...'], ref_s(iter));
            break
        end
        
        time_opt_settletime_s(iter) = U.Time(end) + Nd*sys_sim.Ts;

        traj_iter = struct('X', X, 'U', U);
        opt_trajs_save{iter} = traj_iter;
        if verbose
            change_current_figure(Fig); 
            hold on
            hands(iter) = plot(ref_s(iter), time_opt_settletime_s(iter)*1000, 'xk'); 
        end
    end

    if verbose
        h = plot(ref_s, time_opt_settletime_s*1000, 'LineWidth', 2); 
    end
    
    % If the bisection failed, truncate the results, so that the
    % simulation will run again.
    if status
        step_data.params.ref_s = ref_s(1:iter-1);
        time_opt_settletime_s = time_opt_settletime_s(1:iter-1);
        opt_trajs_save = opt_trajs_save(1:iter-1);
    end
    
    
    results.time_opt_settletime_s =  time_opt_settletime_s;
    results.timeopt_traj_s = opt_trajs_save;
        
    step_data.results = results;
        
    if step_data.savedata
        save(step_data.file, 'step_data')
    end
    
    
end % end main function



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



