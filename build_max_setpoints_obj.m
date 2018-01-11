% max_ref_data = build_max_setpoints(data)

% I want to investigate methods to increase the condition number of the
% Hessian

% Build up the correct model from what is saved from sysID. Ie, put the
% thing in a 

function step_data = build_max_setpoints_obj(step_data, varargin)


% ------------------ Linear + delU Saturation Case ---------------------- %
% Iterate over a bunch of gammas and try to find the maximum setpoint for
% each one. This code is pretty niave. We start at a very low setpoint and
% slowly increase the setpoint until the settling time is reported as NaN.
% A bisection search would probably be faster. This is complicated though
% by the fact that instability seems to occur in the "middle": if the
% setpoint is large enough, we dont have the stability problem, which is
% weird.

% pull data out of data variable:
    


    defaultForce = 0;
    p = inputParser;
    p.addParameter('force', defaultForce);
    p.addParameter('fid', 1);
    parse(p, varargin{:});
    force = p.Results.force;
    fid = p.Results.fid;
    if stepdata_struct_unchanged(step_data) && ~force
        load(step_data.file)
        fprintf(fid, 'LOG (build_max_setpoints, mpc_on=%d)\n', step_data.params.sim_struct.mpc_on);
        fprintf(fid, ['Data appears to be the same. Loading data ',...
                 'without re-calculation.\n\n']);
        return
    end
    fprintf(fid, 'LOG (build_max_setpoints, mpc_on=%d)\n', step_data.params.sim_struct.mpc_on);
    fprintf(fid, 'data has changed: re-building max setpoints.\n');
    
    params = step_data.params;
    % Expose parameters:
    verbose = step_data.verbose;
     
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
            PB = ProgressBar(length(params.gam_s), 'start_str', 'Lin');
        end

        max_setpoints = build_max_sp_local(params, Figs, ...
                                           PB, verbose);
    else
        % Pre-allocate.
        max_setpoints = repmat(0*params.gam_s, length(N_mpc_s), 1);
        sim_struct.mpc_on = true;
        for mpc_iter = 1:length(params.N_mpc_s)
            N_mpc_iter = params.N_mpc_s(mpc_iter);
            if verbose > 0
                start_str = sprintf('MPC, N=%.0f', N_mpc_iter);
                PB = ProgressBar(length(params.gam_s), 'start_str', start_str);
            end

            max_setpoints_iter = build_max_sp_local(params, Figs, PB, verbose, ...
                                                    N_mpc_iter);
            max_setpoints(mpc_iter, :) = max_setpoints_iter;
        end
    
    end

    results.(max_sp_fieldname) = max_setpoints;
    step_data.results = results;
        
    if verbose > 1
        hlin.DisplayName = 'LQR lin + sat';
        xlabel('$\gamma$', 'interpreter', 'latex', 'FontSize', 16)
        ylabel('Max Ref', 'interpreter', 'latex', 'FontSize', 16)
        ylim([0, 5.1])
        title('with delay')
        % hmpc.DisplayName = 'MPC';
        % legend([hlin, hmpc])
    end

    if step_data.savedata
        save(step_data.file, 'step_data')
        if verbose >1
            saveas(F1, step_data.fig_files(1))
            saveas(F1, step_data.fig_files(2))
        end
    end

end % END MAIN FUNCTION

function max_setpoints_iter = build_max_sp_local(params, Figs, PB, verbose, ...
                                                 varargin)
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
        colrs = get(gca, 'colororder');
        PB.upd(0);
    end
    for gam_iter = 1:length(gam_s)
        gamma = gam_s(gam_iter);
        % Update either K or the mpcProb.
        sim_struct = update_sim_struct(sim_struct, params, gamma, N_mpc_iter);
        
        rmax_ind = find(ref_s == ref_max, 1);
        rmax_ind = max([rmax_ind-5, 1]);
        
        fig_base = 10*gam_iter;
        [ref_max, t_settle_s] = find_ref_max(sim_struct, ref_s(rmax_ind:end),...
                                             'verbose', verbose,...
                                             'fig_base', fig_base);
        max_setpoints_iter(gam_iter) = ref_max;
        
        if verbose > 1
            if mod(gam_iter, 10)==0
                fprintf('\ndeleting current figures\n\n')
                for i=1:gam_iter
                    close(figure(10*i + 1))
                    close(figure(10*i + 2))
                    close(figure(10*i + 3))
                end
            end
            change_current_figure(Figs(2));
            k_max = find(t_settle_s == 0, 1, 'first') -1;
            plot(ref_s(1:k_max), t_settle_s(1:k_max)*1e3)
            ylim([0, 10])

            if exist('hlin', 'var')
                delete(hlin)
                drawnow()
            end

            change_current_figure(Figs(1))
            hlin = plot(gam_s(1:gam_iter), max_setpoints_iter(1:gam_iter),...
                        '-o', 'Color', colrs(1,:), 'LineWidth', 2);
            drawnow()
            hold on
        end
        if verbose > 0
            PB.upd(gam_iter);
        end
        
        %     keyboard
    end % end MAIN LOOP 
    
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




% function status = data_struct_unchanged(data_struct)
% % This function will return 1 if
%     % (a) data_struct.file exists
%     % (b) the fields .Q, .sys, .ref_s, .gam_s are the same
%     %     as the data contained in data_struct.file
%     % (c) the field (.max_setpoints_lin)(.max_setpoints_mpc) exists,
%     %     depending on the flag data_struct.mpc_on
    
%     if ~exist(data_struct.file, 'file')
%         status = 0;
%         % keyboard
%         return 
%     else
%         load(data_struct.file)
%     end
%     fields1 = {'Q', 'ref_s', 'gam_s'};
%     fields2 = {'A', 'B', 'C', 'D', 'Ts'};
%     if fields_have_changed(data_struct, max_ref_data, fields1)
%         status = 0;
%         % keyboard
%         return
%     elseif isfield(data_struct, 'sys') && isfield(max_ref_data, 'sys')
%         if fields_have_changed(data_struct.sys, max_ref_data.sys, fields2)
%             status = 0;
%             % keyboard
%         end
%     else
%         % only way we get here is if .sys is not afield of either structure.
%         status = 0;
%         keyboard
%         return
%     end
            
%     if data_struct.mpc_on == 0
%         if ~isfield(max_ref_data, 'max_setpoints_lin')
%             status = 0;
%             % keyboard
%             return
%         end
%     elseif data_struct.mpc_on == 1
%         if ~isfield(max_ref_data, 'max_setpoints_mpc')
%             status = 0;
%             % keyboard
%             return
%         end
%         if ~isequal(data_struct.N_mpc_s, max_ref_data.N_mpc_s)
%             status = 0;
%             % keyboard
%             return
%         end
%     end
%     status = 1;
% end


% function status = fields_have_changed(strct1, strct2, fields)
%     for fld = fields
%         A = isfield(strct2, fld{1}) || isprop(strct2, fld{1} );
%         B = isfield(strct1, fld{1}) || isprop(strct1, fld{1} );
%         if A && B
%             if ~isequal(strct2.(fld{1}), strct1.(fld{1}))
%                 % keyboard
%                 status = 1;
%                 return
%             end
%         else
%             status = 1;
%             % keyboard
%             return 
%         end
%     end
%     status = 0;
% end
