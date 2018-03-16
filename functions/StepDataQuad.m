classdef StepDataQuad < StepData
    
    properties
        % params;
        % file;
        % fig_files;
        % savedata;
        % results;
        max_ref_judge; % Handle to a function which decides if the
                       % maximum reference has been found. This
                       % will be passed into find_max_ref;
        ts_by_rmax_results;
    end
    
    methods
        function self = StepDataQuad(Params, varargin)
            self = self@StepData(Params, varargin{:})
            self.max_ref_judge = @max_sp_judge_default;
            ts_by_rmax_results = {};
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
            
            ref_max_idx_prior = 1;
            max_setpoints_idx = 0*gam_s;
            max_recommended_sps_idx = 0*gam_s;
            PB.upd(0);
            
            data = cell(1, length(gam_s));
            for gam_iter = 1:length(gam_s)
                gamma = gam_s(gam_iter);
                % Update either K or the mpcProb.
                sim_struct = self.update_sim_struct(sim_struct, gamma);

                % Warm starting
                % rmax_ind = find(ref_s == ref_max_prior, 1);
                % rmax_ind = max([rmax_ind-5, 1]);
                rmax_ind = 1;

                max_ref_judge = self.max_ref_judge.build_ref_max_judge_dynamic(gamma); %#ok<PROPLC>
                
                fig_base = 10*gam_iter;
                data_iter = find_ref_max(sim_struct, ref_s(rmax_ind:end),...
                                         'verbose', verbose,...
                                         'fig_base', fig_base,...
                                         'max_sp_judge', max_ref_judge); %#ok<PROPLC>
                
                max_setpoints_idx(gam_iter) = data_iter.ref_max_idx;
                max_recommended_sps_idx(gam_iter) = data_iter.ref_max_recommended_idx;
                
                % ref_max_idx_prior = data_iter.ref_max_idx;
                data{gam_iter} = data_iter;
                
                % keyboard
                if verbose > 1
                    plot_local(data_iter.t_settle_s, gam_s, ref_s, max_recommended_sps_idx,...
                               max_setpoints_idx, gam_iter, self.logger);
                end
                PB.upd(gam_iter);
                
                %     keyboard
            end % end MAIN LOOP 
                result_s.max_setpoints_idx = max_setpoints_idx;
                result_s.max_setpoints = ref_s(max_setpoints_idx);
                result_s.max_recommended_sps_idx = max_recommended_sps_idx;
                result_s.max_recommended_sps = ref_s(max_recommended_sps_idx);
                result_s.data = data;
        end
        
        function [ax, h] = plot_ts_by_ref_max_judge(self, rmax, ...
                                                          exp_idx, varargin)
            if isempty(self.ts_by_rmax_results)
                error(['You must generate the ts by rmax results ' ...
                       'first. Run self. ts_by_ref_max(ref_max, exp_idx'])
            end
            if isempty(exp_idx)
                exp_idx = 1;
            end
            if isempty(varargin)
                f = figure;
                ax = gca;
            else
                ax = varargin{1};
                varargin(1) = [];
            end
            
            h = plot(ax, self.ts_by_rmax_results{exp_idx}.ref_s, ...
                 self.ts_by_rmax_results{exp_idx}.t_settle_s*1e3,...
                     varargin{:});
            h.DisplayName = sprintf('exp-idx = %.0f, rmax-des = %.2f, $\\gamma=%.0f$', ...
                                    exp_idx,...
                                    self.ts_by_rmax_results{exp_idx}.rmax_des,...
                                    self ...
                                    .ts_by_rmax_results{exp_idx}.gamma);
            xlabel('setpoint')
            ylabel('settle time [ms]')
        end
        
        function [ax, h] = plot_ts_perc_increase_by_rmax(self,rmax, ...
                                                         exp_idx, ts_other, varargin)
        % [ax, h] = plot_ts_perc_increase_by_rmax(self,rmax, exp_idx, ts_other, varargin)
        % 
        % Plot percent increase for a maximum desired reference
        % compared to data in ts_other
        % 
        % Inputs
        % -----
        %  rmax : 
        %  
        %  exp_idx : for linear, this is 1. For MPC, this is the
        %  N_mpc index in self.params.N_mpc_s
        %  ts_other : a vector of settling times to compare to.
        %
            
            if isempty(self.ts_by_rmax_results)
                error(['You must generate the ts by rmax results ' ...
                       'first. Run self. ts_by_ref_max(ref_max, exp_idx'])
            end
            if isempty(exp_idx)
                exp_idx = 1;
            end
            if isempty(varargin)
                f = figure;
                ax = gca;
            else
                ax = varargin{1};
                varargin(1) = [];
            end
            ts_self = self.ts_by_rmax_results{exp_idx}.t_settle_s;
            ref_s = self.ts_by_rmax_results{exp_idx}.ref_s;
            
            len_ts_self = length(ts_self);
            
            if len_ts_self > length(ts_other)
                error(['The number of settle times in ' ...
                       'self.ts_by_rmax_results must be at least ' ...
                       'length of ts_other']);
            end
            
            ts_other = ts_other(1:len_ts_self);
            perc_increase = (ts_self./ts_other)*100;
            h = plot(ax, ref_s, perc_increase, varargin{:});
            h.DisplayName = sprintf('exp-idx = %.0f, rmax-des = %.2f, $\\gamma=%.0f$', ...
                                    exp_idx,...
                                    self.ts_by_rmax_results{exp_idx}.rmax_des,...
                                    self.ts_by_rmax_results{exp_idx}.gamma);
            xlabel('setpoint')
            ylabel('settle time \% increase')
        end
        
        function self = ts_by_ref_max(self, rmax, exp_idx)
        % self = ts_by_ref_max(self, rmax, exp_idx)    
            if isempty(self.results)
               error(['You must first load or generate results before  ' ...
                      'this function can be used'])
            end
            
            ref_s  = self.params.ref_s;
            gam_s =  self.params.gam_s;
            max_recommended_sps = self.results{exp_idx}.max_recommended_sps;
            max_recommended_sps_idx = self.results{exp_idx}.max_recommended_sps_idx;

            gam_idx = find(max_recommended_sps >=rmax, 1, 'first');

            result.gamma = gam_s(gam_idx);
            ref_max_idx = max_recommended_sps_idx(gam_idx);
            
            t_settle_s = self.results{exp_idx}.data{gam_idx}.t_settle_s;
            
            result.ref_s = ref_s(1:ref_max_idx);

            % we typically generate ALL of the trajs, so have ALL t_settle_s
            result.t_settle_s = t_settle_s(1:ref_max_idx);
            result.rmax_des = rmax;
            self.ts_by_rmax_results{exp_idx} = result;
        end
        
        function sim_struct = update_sim_struct(self, sim_struct, gamma)
        % sim_struct = update_sim_struct(self, sim_struct, gamma)
        % 
        % Build sim_struct based off the current data in sim_struct and
        % control weight gamma. If sim_struct.N_mpc ==0, build the feedback
        % gain K_lqr. Otherwise, build an condensedMPCprob_OA
        % 
            warning('off', 'MATLAB:structOnObject');
            params = self.params;
            if sim_struct.N_mpc == 0 % we are in linear mode then.
                sim_struct.K_lqr = dlqr(params.sys.a, params.sys.b,...
                                        params.Q, gamma, self.params.S);
            else % if N_MPC_iter not empty --> in MPC mode.
                du_max = sim_struct.du_max;
                N_mpc  = sim_struct.N_mpc;
                Qp = dare(params.sys.a, params.sys.b, params.Q, gamma,...
                    self.params.S); 
                sim_struct.mpcProb1 = condensedMPCprob_OA(params.sys, N_mpc,...
                                                       params.Q, Qp, gamma,...
                                                       self.params.S);
                CON = CondenCon([], [], N_mpc);
                CON.add_input_con('box', [-du_max, du_max]);
                sim_struct.mpcProb1.CON = CON;
                %sim_struct.mpcProb1.add_U_constraint('box', [-du_max, du_max]);
                sim_struct.K_lqr = params.sys.C*0;
            end
        end
        
    end
    
end


function plot_local(t_settle_s, gam_s, ref_s, max_recommended_sps_idx,...
                    max_setpoints_idx, gam_iter, logger)
% expose for easy access while plotting:
    
    if mod(gam_iter, 10)==0
        logger('\ndeleting current figures\n\n');
        for i=1:gam_iter
            close(figure(10*i + 1));
            close(figure(10*i + 2));
            close(figure(10*i + 3));
        end
    end
    max_setpoints = ref_s(max_setpoints_idx(1:gam_iter));
    max_setpoints_recommended = ref_s(max_recommended_sps_idx(1:gam_iter));
    F1 = figure(1); hold on;
    xlabel('$\gamma$', 'interpreter', 'latex', 'FontSize', 16)
    ylabel('Max Ref', 'interpreter', 'latex', 'FontSize', 16)
    
    F2 = figure(2);  hold on;
    ylabel('settle time [ms]')
    xlabel('ref_f')

    change_current_figure(F2);
    % colrs = get(gca, 'colororder');
    k_max = find(t_settle_s ~= 0, 1, 'last');
    % keyboard;
    k_max_rec = max_recommended_sps_idx(gam_iter);
    
    plot(ref_s(1:k_max), t_settle_s(1:k_max)*1e3)
    plot(ref_s(k_max_rec), t_settle_s(k_max_rec)*1e3, 'xk')
    % ylim([0, 10])
    hold on
    drawnow()

    change_current_figure(F1)
    hlin = plot(gam_s(1:gam_iter), max_setpoints(1:gam_iter),...
                '-o',  'LineWidth', 1);
    plot(gam_s(1:gam_iter), max_setpoints_recommended(1:gam_iter), '-x',...
         'LineWidth', 2);

    drawnow()
    hold on
end


function ref_max_recommended_found = max_sp_judge_default(varargin)
    ref_max_recommended_found = false;
end

