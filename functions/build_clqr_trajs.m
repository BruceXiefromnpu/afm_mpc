% clqr_data = build_clqr_trajs(step_data, varargin)
%
% Builds CLQR trajectories parameteriezed by gamma. Will run a sequence
% of simulations over a list of references, for each gamma. Will record
% both the resulting trajectories as well as the associated settling
% times. 
%
% Inputs
% -----
%   step_data : a class instance of StepData 
% 
% Outputs
%  -----
%  step_data : the same class that was passed in, but with the
%  results field populated. The results field shall have the form:
%  
%    N.B. The results fields are cell arrays, to account for the
%    occasion that we gam_s is a vector.
% 
%    results.settle_times_opt_cell : settle times for each reference.
% 
%    results.opt_trajs_cell        : the trajectories associated
%    with each reference (Y, U, X).
%
%
% The big function of this function is to check if step_data.file
% exists. If it does, then the function compares the data saved in that
% .mat file to the data provided in data_struct. If the data is the same,
% the function simply loads the .mat file and returns the structure
% contained therin. If the data has changed, the function will perform
% the simulations again, saving the new data into the .mat file specified
% in data_struct.file.
%
%
% See Also: opt_traj_gen, StepDataCLQR

function step_data = build_clqr_trajs(step_data, varargin)
% Description.
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
        fprintf(fid, 'LOG: (build_clqr_trajs)\n');
        fprintf(fid, ['Data appears to be the same. Loading data ',...
                 'without re-calculation.\n\n']);
        return
    end
    fprintf(fid, 'LOG (build_clqr_trajs\n');
    fprintf(fid, 'data has changed: re-building max setpoints.\n');
    
    % keyboard
    params = step_data.params;
    
    % Expose parameters:
    gam_s = params.gam_s;
    ref_s = params.ref_s;
    du_max = params.du_max;
    N_traj = params.N_traj;
    mpc_mode = params.mpc_mode;
   
    % Pre-allocate
    settle_times_opt_save = cell(1, length(gam_s));
    opt_trajs_save = cell(1, length(gam_s));

    
    if verbose
        Fig = figure(200);
        hands = [];
        change_current_figure(Fig);
        colrs = get(gca, 'colororder');
    end
    for iter = 1:length(gam_s)
        gamma = gam_s(iter);

        [traj_s, settle_times_opt] = opt_traj_gen(params.Q, gamma, N_traj,...
                                                  params.sys, ...
                                                  ref_s, 0, 'uMax', du_max, 'verbose', 0,...
                                                  'mpc_mode', mpc_mode);
        settle_times_opt_save{iter} = settle_times_opt;
        opt_trajs_save{iter} = traj_s;
        if verbose
            change_current_figure(Fig); 
            hands(iter) = plot(ref_s, settle_times_opt_save{iter}*1000, 'LineWidth', 2); 

            set(hands(iter), 'DisplayName', sprintf('$\\gamma = %.0f$', gamma),...
                             'Color', colrs(iter, :));
            leg = legend(hands);
            set(leg, 'interpreter', 'latex');
            ylabel('settle time [ms]', 'FontSize', 16)
            xlabel('setpoint', 'FontSize', 16)
            drawnow()
            grid on
        end
    end

    
    results.settle_times_opt_cell =  settle_times_opt_save;
    results.opt_trajs_cell = opt_trajs_save;
        
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


