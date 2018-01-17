classdef StepDataQuad_parfor < StepData
    
    properties
        % params;
        % file;
        % fig_files;
        % savedata;
        % results;
    end
    
    methods
        function self = StepDataQuad_parfor(Params, varargin)
            self = self@StepData(Params, varargin{:})
        end
        
        function result_s = build_max_sp_local(self, PB, verbose, varargin)
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

        % Initialize the "warm starting"
            sim_struct = self.params.sim_struct;
            if ~isempty(varargin)
                sim_struct.N_mpc = varargin{1};
            end
            gam_s = self.params.gam_s;
            ref_s = self.params.ref_s;
            
            ref_max_prior = ref_s(1);
            max_setpoints_iter = 0*gam_s;
            PB.upd(0);
            
            data = cell(1, length(gam_s));
            for gam_iter = 1:length(gam_s)
                gamma = gam_s(gam_iter);
                % Update either K or the mpcProb.
                sim_struct = self.update_sim_struct(sim_struct, gamma);
                % keyboard 
                % Warm starting
                rmax_ind = find(ref_s == ref_max_prior, 1);
                rmax_ind = max([rmax_ind-5, 1]);
                
                fig_base = 10*gam_iter;
                data_iter = find_ref_max(sim_struct, ref_s(rmax_ind:end),...
                                         'verbose', verbose,...
                                         'fig_base', fig_base);
                max_setpoints(gam_iter) = data_iter.ref_max;
                ref_max_prior = data_iter.ref_max;
                data{gam_iter} = data_iter;
                
                if verbose > 1
                    plot_local(data_iter.t_settle_s, gam_s, ref_s, max_setpoints, ...
                               gam_iter, self.logger);
                end
                PB.upd(gam_iter);
                
                %     keyboard
            end % end MAIN LOOP 
                result_s.max_setpoints = max_setpoints;
                result_s.data = data;
        end

        function sim_struct = update_sim_struct(self, sim_struct, gamma)
            params = self.params;
            if sim_struct.N_mpc == 0 % we are in linear mode then.
                sim_struct.K_lqr = dlqr(params.sys.a, params.sys.b,...
                                        params.Q, gamma);
            else % if N_MPC_iter not empty --> in MPC mode.
                du_max = sim_struct.du_max;
                N_mpc  = sim_struct.N_mpc;
                Qp = dare(params.sys.a, params.sys.b, params.Q, gamma); 
                sim_struct.mpcProb1 = condensedMPCprob_OA(params.sys, N_mpc,...
                                                       params.Q, Qp, gamma);
                sim_struct.mpcProb1.add_U_constraint('box', [-du_max, du_max]);
                sim_struct.K_lqr = params.sys.C*0;
            end
        end
        
    end
    
end


function plot_local(t_settle_s, gam_s, ref_s, max_setpoints, gam_iter, ...
                    logger)
% expose for easy access while plotting:
    
    if mod(gam_iter, 10)==0
        logger('\ndeleting current figures\n\n');
        for i=1:gam_iter
            close(figure(10*i + 1));
            close(figure(10*i + 2));
            close(figure(10*i + 3));
        end
    end
    F1 = figure(1); hold on;
    xlabel('$\gamma$', 'interpreter', 'latex', 'FontSize', 16)
    ylabel('Max Ref', 'interpreter', 'latex', 'FontSize', 16)
    
    F2 = figure(2);  hold on;
    ylabel('settle time [ms]')
    xlabel('ref_f')

    change_current_figure(F2);
    % colrs = get(gca, 'colororder');
    k_max = find(t_settle_s ~= 0, 1, 'last');
    
    plot(ref_s(1:k_max), t_settle_s(1:k_max)*1e3)
    % ylim([0, 10])
    hold on
    drawnow()

    change_current_figure(F1)
    hlin = plot(gam_s(1:gam_iter), max_setpoints(1:gam_iter),...
                '-o',  'LineWidth', 2);

    drawnow()
    hold on
end
