classdef MaxSpJudgeCLQR

    properties
        ref_s;
        ts_s_compare;
        thresh_hold;
        step_data_clqr;
    end
    
    methods
        function self = MaxSpJudgeCLQR(step_data_clqr, thresh_hold)
            if ~isa(step_data_clqr, 'StepDataCLQR')
                error(['first parameter must be of class ' ...
                       'StepDataCLQR'])
            end
            
            self.ref_s = step_data_clqr.params.ref_s;
            self.ts_s_compare = ...
                step_data_clqr.results.settle_times_opt_cell{1};
            self.thresh_hold = thresh_hold;
            
            step_data_clqr.results = []; % Dont carry all that
                                         % data. 
            self.step_data_clqr = step_data_clqr;

        end
        
        function ref_max_judge = build_ref_max_judge_single_baseline(self, gam_current)
                ref_max_judge = @self.max_sp_judge_single_baseline;
        end
            
        function ref_max_judge = build_ref_max_judge_dynamic(self, ...
                                                             gam_current)
        % This is a factory function that builds a function handle to the
        % function max_sp_judge_dynamic based on the current gamma. 
            self.step_data_clqr.params.gam_s = gam_current;
            
            % Expose parameters and build CLQR problem for current gamma:
            
            params = self.step_data_clqr.params;
            sys = params.sys;
            du_max = params.du_max;
            N_traj = params.N_traj;
            mpc_mode = params.mpc_mode;
            Q = params.Q;
            S = params.S;
            Qp = dare(sys.a, sys.b, Q, gam_current, S);
            
            NLQR_prob = condensedMPCprob_OA(sys, N_traj, Q, Qp, gam_current, S);

            CON = CondenCon([], [], NLQR_prob.N_mpc);
            CON.add_input_con('box', du_max);
            NLQR_prob.CON = CON;
            
            % NLQR_prob.add_U_constraint('box', du_max);
            Nx = SSTools.getNxNu(sys);
            
            % Build the function handle:
            ref_max_judge = @(ref_f, idx, ts_other, Y)self.max_sp_judge_dynamic(ref_f, idx, ...
                             ts_other, Y, NLQR_prob, Nx, sys);
        end
        function ref_max_recommended_found = max_sp_judge_dynamic(self,...
            ref_f, idx, ts_other, Y, NLQR_prob, Nx, sys)
            % Assume zero initial condition.
            x0_err = 0*Nx - ref_f*Nx;
            [U, Xerr] = NLQR_prob.solve(x0_err, 'getX', true); 
            X = Xerr + ref_f*Nx;
            Y = sys.C*X;
            tvec = (0:1:size(U, 2)-1)*sys.Ts; 
            ts_self = settle_time(tvec, Y, ref_f, 0.01*ref_f, [], [], 30);
            perc_increase = (ts_other/ts_self)*100;
            
            if perc_increase > self.thresh_hold || isnan(perc_increase)
                ref_max_recommended_found = true;
            else
                ref_max_recommended_found = false;
            end
            
        end
        
        function ref_max_recommended_found = max_sp_judge_single_baseline(self,...
            ref_f, idx, ts_other, Y)
            % keyboard
            if self.ref_s(idx) ~= ref_f
                error(['It seems something is amiss: ref_f must be ' ...
                       'the same to compare.']);
            end
            ts_self = self.ts_s_compare(idx);
            
            perc_increase = (ts_other/ts_self)*100;
            
            if perc_increase > self.thresh_hold || isnan(perc_increase)
                ref_max_recommended_found = true;
            else
                ref_max_recommended_found = false;
            end
        end
    end  % methods

end % classdef