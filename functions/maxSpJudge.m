classdef maxSpJudge

    properties
        ref_s;
        ts_s_compare;
        thresh_hold;
    end
    
    methods
        
        function self = MaxSpJudge(step_data, thresh_hold)
            self.ref_s = step_data.params.ref_s;
            self.ts_s_compare = ...
                step_data.results.settle_times_opt_cell{1};
            self.thresh_hold = threshold;
        end
        
        function ref_max_recommended_found = max_sp_judge(self, ref_f, idx, ts_other, Y)
            
            if self.ref_s(iter) ~= ref_f
                error(['It seems something is amiss: ref_f must be ' ...
                       'the same to compare.']);
            end
            ts_self = self.ts_s_compare(iter);
            
            perc_increase = (ts_other/ts_self)*100;
            
            if perc_increase > self.thresh_hold || isnan(perc_increase)
                ref_max_recommended_found = true;
            else
                ref_max_recommended_found = false;
            end
            
            
        end
    end  % methods

end % classdef