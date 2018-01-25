classdef MaxSpJudge

    properties
        ref_s;
        ts_s_compare;
        thresh_hold;
    end
    
    methods
        function self = MaxSpJudge(step_data_clqr, thresh_hold)
            if ~isa(step_data_clqr, 'StepDataCLQR')
                error(['first parameter must be of class ' ...
                       'StepDataCLQR'])
            end
            
            self.ref_s = step_data_clqr.params.ref_s;
            self.ts_s_compare = ...
                step_data_clqr.results.settle_times_opt_cell{1};
            self.thresh_hold = thresh_hold;
        end
        
        function ref_max_recommended_found = max_sp_judge(self, ref_f, idx, ts_other, Y)
            
            if self.ref_s(idx) ~= ref_f
                error(['It seems something is amiss: ref_f must be ' ...
                       'the same to compare.']);
            end
            ts_self = self.ts_s_compare(idx);
            
            perc_increase = (ts_other/ts_self)*100;

            
            % fprintf('perc increase = %f\n', perc_increase);
            
            if perc_increase > self.thresh_hold || isnan(perc_increase)
                ref_max_recommended_found = true;
            else
                ref_max_recommended_found = false;
            end
            
            
        end
    end  % methods

end % classdef