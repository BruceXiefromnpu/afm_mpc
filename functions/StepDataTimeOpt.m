classdef StepDataTimeOpt < StepData
    
    properties
    end
    
    methods
        function self = StepDataTimeOpt(Params, varargin)
            self = self@StepData(Params, varargin{:})
            % p = inputParser;
            % p.addParameter('savedata', true)
            % p.addParameter('file', '')
            % p.addParameter('fig_files', '')
            % parse(p, varargin{:});
            
            
            % self.params = Params;
            % self.file = p.Results.file;
            % self.fig_files = p.Results.fig_files;
            % self.savedata =  p.Results.savedata;
            % self.results = [];
        end
        
        function plot_single_traj(self, index, ax)
            if ~exist('ax', 'var')
                ax = gca();
            end
            ref_s = self.params.ref_s;
            clqr_settletime_s = self.results.settle_times_opt_cell{1};
            h = plot(ax, ref_s, clqr_settletime_s*1000, varargin{:})
            set(h, 'DisplayName', 'CLQR')
        end
        
    end
    
end
