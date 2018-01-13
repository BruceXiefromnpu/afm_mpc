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
        
        function plot_ref_vs_settle(self,ax, varargin)
        % plot_ref_vs_settle(self,ax, varargin)
        % plot reference vs settle time for data contained in
        % self.data. 
        %   -- If ax is empty, will plot to gca().
        %   -- varargin is passed straight to matlabs plot function. 
            if ~exist('ax', 'var')
                ax = gca();
            elseif isempty(ax)
                ax = gca();
            end
            ref_s = self.params.ref_s;
            clqr_settletime_s = self.results.settle_times_opt_cell{1};
            h = plot(ax, ref_s, clqr_settletime_s*1000, varargin{:});
            set(h, 'DisplayName', 'CLQR')
        end
        
        function [h, ax] = plot_single_ytraj(self, index, ax, ...
                                             varargin)
        % plot the y-trajectory held at self.results.opt_trajs_cell{1}.Y_vec_s{index}
        %   -- If ax is empty, will plot to gca().
        %   -- varargin is passed straight to matlabs plot function. 

            if ~exist('ax', 'var')
                ax = gca();
            elseif isempty(ax)
                ax = gca();
            end
            traj_y = self.results.opt_trajs_cell{1}.Y_vec_s{index}
            hy = plot(ax, traj_y.Time, traj_y.Data, varargin{:});
            
        end

        function [h, ax] = plot_single_utraj(self, index, ax, varargin)
        % Plot the u-trajectory held at self.results.opt_trajs_cell{1}.U_vec_s{index}
        %   -- If ax is empty, will plot to gca().
        %   -- varargin is passed straight to matlabs plot function. 
            if ~exist('ax', 'var')
                ax = gca();
            elseif isempty(ax)
                ax = gca();
            end

            traj_u = self.results.opt_trajs_cell{1}.U_vec_s{index}
            hu = plot(ax, traj_u.Time, traj_u.Data, varargin{:});
            
        end
        
    end
    
end
