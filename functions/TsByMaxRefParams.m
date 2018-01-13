classdef TsByMaxRefParams 
    properties
        StepData;
        rmax_s;
        file;
        results;
    end
    methods
        function self = TsByMaxRefParams(StepData, rmax_s)
            self.StepData = StepData;
            self.rmax_s = rmax_s;
        end
        function self = run_ts_by_refs(self)
            
            results = cell(1, length(self.rmax_s));
            for rmax_iter = 1:length(self.rmax_s)
                rmax = self.rmax_s(rmax_iter);
                % Expose for easy access.
                ref_s = self.StepData.params.ref_s;
                gam_s = self.StepData.params.gam_s;
                sim_struct = self.StepData.params.sim_struct;
                max_setpoints_by_gamma = ...
                    self.StepData.results.max_setpoints_lin;
                % For this size reference, find the first gamma we
                % can use. 
                kk = find(max_setpoints_by_gamma >=rmax, 1, 'first');
                gamma = gam_s(kk);
                if rmax == 4.8
                    gamma = 9035;
                end
                
                
                
                % Get the feedback gain that can do it.
                K_lqr = dlqr(self.StepData.params.sys.a, self.StepData.params.sys.b,...
                             self.StepData.params.Q, gamma);
                
                
                
                % Index to slice ref_s from the first, up the one
                % corresponding to our max setpoints. N.B.
                % max_setpoints and ref_s ARE NOT equal!.
                k_ref = find(ref_s == max_setpoints_by_gamma(kk), 1, ...
                             'first');
                t_settle_s = zeros(1, k_ref);
                % for each ref BELOW rmax, simulate and compute the
                % settling times:
                
                for k=1:k_ref
                    mpc_on = 0;
                    ref_f = ref_s(k);
                    y1 = sim_local(sim_struct, K_lqr, ref_f);
                    
                    [t_settle, k_s] = settle_time(y1.time, y1.Data, ref_f, 0.01*ref_f,...
                                                           [], [], 30);
                    t_settle_s(k) = t_settle;
                end
                result_rmx_iter.t_settle_s = t_settle_s;
                result_rmx_iter.ref_s_to_rmax = ref_s(1:k);
                result_rmx_iter.gamma = gamma;
                results{rmax_iter} = result_rmx_iter;
            end
            
            self.results = results;
        end
        
        function [ax, hands, leg] = plot_ts_v_r2max(self, ax)
            if isempty(self.results)
                error('No Results to plot!')
            end
            
            if isempty(ax)
                ax = gca();
            end
            hold on
            hands = [];
            % keyboard
            for k=1:length(self.results)
                rmax = self.results{k}.ref_s_to_rmax(end);
                gamma = self.results{k}.gamma;
                leg_str = sprintf('$\\gamma = %.0f$, $r_{max} = %.1f$', gamma, rmax);
                h = plot(ax, self.results{k}.ref_s_to_rmax, ...
                         self.results{k}.t_settle_s*1000);
                set(h, 'DisplayName', leg_str);
                xlabel('reference', 'FontSize', 14, 'interpreter', ...
                       'latex')
                ylabel('settle time [ms]', 'FontSize', 14, 'interpreter', ...
                       'latex')
                hands(k) = h;
            end
            leg = legend(hands);
            set(leg, 'interpreter', 'latex', 'FontSize', 14, ...
                     'Location', 'NorthWest');
            
            
        end
        
        function save(self)
            save(self.file, 'self')
        end
        function eqq = isequal(self, other, varargin)
            if length(varargin > 0)
                strict = 1;
            else
                strict = 0;
            end
                        
            if strict
                eqq = isequal(self, other);
            else
                eqq1 = isequal(self.StepData.params, ...
                               other.StepData.params);
                eqq2 = isequal(self.rmax_s, other.rmax_s)
                eqq = eqq1 && eqq2;
            end
        end
    end
end

function [Y, U, dU] = sim_local(sim_struct, K_lqr, ref_f)
    options = simset('SrcWorkspace','current');
    % Expose the sime struct to simulink.
    PLANT = sim_struct.PLANT;
    trun = sim_struct.trun;
    mpcProb1 = sim_struct.mpcProb1;
    du_max = sim_struct.du_max;
    mpc_on = sim_struct.mpc_on;
    xss = sim_struct.xss;
    Nx = sim_struct.Nx;
    
    x0 = sim_struct.xss*0;
    uss_0 = 0;
    Ts = PLANT.Ts;
    sim('MPC_fp', [], options)
    
    Y = y_mpcDist;
    U = u_mpcDist;
    dU = du_mpcDist;
end
